#include "json_value.hpp"
#include "json.hpp"
#include "parse.hpp"
#include "types.hpp"
#include <cctype>
#include <cmath>
#include <iomanip>

namespace kutta {
namespace {
class Parser {
  std::string_view text;std::size_t at=0;
  [[noreturn]] void fail(std::string_view why) const {throw Error("JSON byte "+std::to_string(at)+": "+std::string(why));}
  char peek() const {return at<text.size()?text[at]:'\0';}
  void space(){while(peek()==' ' || peek()=='\n' || peek()=='\r' || peek()=='\t')++at;}
  void expect(char c){if(peek()!=c)fail("unexpected character");++at;}
  unsigned hex4(){unsigned n=0;for(int i=0;i<4;++i){const char c=peek();++at;n*=16;if(c>='0' && c<='9')n+=c-'0';else if(c>='a' && c<='f')n+=c-'a'+10;else if(c>='A' && c<='F')n+=c-'A'+10;else fail("invalid Unicode escape");}return n;}
  void utf8(std::string& s,unsigned n){if(n<128)s+=char(n);else if(n<2048){s+=char(0xc0|(n>>6));s+=char(0x80|(n&63));}else if(n<65536){s+=char(0xe0|(n>>12));s+=char(0x80|((n>>6)&63));s+=char(0x80|(n&63));}else{s+=char(0xf0|(n>>18));s+=char(0x80|((n>>12)&63));s+=char(0x80|((n>>6)&63));s+=char(0x80|(n&63));}}
  std::string string(){
    expect('"');std::string out;
    while(peek()!='"') {
      const auto c=static_cast<unsigned char>(peek());if(c<32)fail("control character or unterminated string");++at;
      if(c!='\\'){out+=char(c);continue;}
      const char e=peek();++at;
      switch(e){case '"':out+='"';break;case '\\':out+='\\';break;case '/':out+='/';break;case 'b':out+='\b';break;case 'f':out+='\f';break;case 'n':out+='\n';break;case 'r':out+='\r';break;case 't':out+='\t';break;
        case 'u':{auto n=hex4();if(n>=0xd800 && n<=0xdbff){expect('\\');expect('u');const auto low=hex4();if(low<0xdc00 || low>0xdfff)fail("invalid surrogate pair");n=0x10000+((n-0xd800)<<10)+(low-0xdc00);}else if(n>=0xdc00 && n<=0xdfff)fail("unpaired surrogate");utf8(out,n);break;}
        default:fail("invalid string escape");}
    }
    ++at;return out;
  }
  Json read(int depth){
    if(depth>64)fail("nesting limit exceeded");space();
    if(peek()=='"')return string();
    if(peek()=='{' || peek()=='['){const bool object=peek()=='{';++at;space();Json::Object obj;Json::Array arr;if(peek()!=(object?'}':']'))for(;;){std::string key;if(object){if(peek()!='"')fail("object key must be a string");key=string();space();expect(':');}auto v=read(depth+1);if(object){if(!obj.emplace(key,std::move(v)).second)fail("duplicate object key: "+key);}else arr.push_back(std::move(v));space();if(peek()!=',')break;++at;space();}expect(object?'}':']');return object?Json(std::move(obj)):Json(std::move(arr));}
    for(auto [token,value]:{std::pair<std::string_view,Json>{"true",true},{"false",false},{"null",nullptr}})if(text.substr(at,token.size())==token){at+=token.size();return value;}
    const auto start=at;if(peek()=='-')++at;
    if(peek()=='0')++at;else {if(peek()<'1' || peek()>'9')fail("expected JSON value");while(peek()>='0' && peek()<='9')++at;}
    if(peek()=='.'){++at;const auto begin=at;while(peek()>='0' && peek()<='9')++at;if(at==begin)fail("fraction needs digits");}
    if(peek()=='e' || peek()=='E'){++at;if(peek()=='+' || peek()=='-')++at;const auto begin=at;while(peek()>='0' && peek()<='9')++at;if(at==begin)fail("exponent needs digits");}
    const auto n=parse_double(text.substr(start,at-start));if(!n || !std::isfinite(*n))fail("number must be finite");return *n;
  }
 public:
  explicit Parser(std::string_view s):text(s){}
  Json parse(){auto v=read(0);space();if(at!=text.size())fail("trailing text");return v;}
};
void dump(std::ostream& s,const Json& j){
  std::visit([&](const auto& v){using T=std::decay_t<decltype(v)>;if constexpr(std::is_same_v<T,std::nullptr_t>)s<<"null";else if constexpr(std::is_same_v<T,bool>)s<<(v?"true":"false");else if constexpr(std::is_same_v<T,double>)json_number(s,v);else if constexpr(std::is_same_v<T,std::string>)s<<json_string(v);else {s<<(std::is_same_v<T,Json::Object>?'{':'[');bool first=true;for(const auto& entry:v){if(!first)s<<',';first=false;if constexpr(std::is_same_v<T,Json::Object>){s<<json_string(entry.first)<<':';dump(s,entry.second);}else dump(s,entry);}s<<(std::is_same_v<T,Json::Object>?'}':']');}},j.value);
}
}
Json parse_json(std::string_view text){KUTTA_REQUIRE(text.size()<=1024*1024,"setup JSON exceeds 1 MiB");return Parser(text).parse();}
std::string dump_json(const Json& j){std::ostringstream s;s.imbue(std::locale::classic());s<<std::setprecision(17);dump(s,j);return s.str();}
}
