#include "kutta/case/sexpr.hpp"
#include <charconv>
#include <cmath>
#include <cctype>
#include <iomanip>
#include <limits>
#include <sstream>

namespace kutta::fluent {
SExpr SExpr::symbol(std::string text) { SExpr x; x.kind=Kind::symbol; x.text=std::move(text); return x; }
std::optional<double> SExpr::number() const {
  if(kind==Kind::integer) return static_cast<double>(integer);
  if(kind==Kind::real) return real;
  return {};
}
std::optional<std::string_view> SExpr::atom_text() const {
  if(kind==Kind::symbol || kind==Kind::string) return text;
  return {};
}
bool SExpr::operator==(const SExpr& b) const {
  return kind==b.kind && text==b.text && integer==b.integer && real==b.real && boolean==b.boolean &&
    elements==b.elements && ((!tail&&!b.tail)||(tail&&b.tail&&*tail==*b.tail));
}
ParseError::ParseError(std::size_t pos, std::string msg)
  :std::runtime_error(std::move(msg)+" at byte "+std::to_string(pos)),offset(pos) {}
namespace {
class Parser {
  std::string_view s;
  ParseLimits limits;
  std::size_t p=0, nodes=0;
  [[noreturn]] void fail(std::string msg) const { throw ParseError(p,std::move(msg)); }
  bool space(char c) const { return std::isspace(static_cast<unsigned char>(c)); }
  bool delimiter(char c) const { return space(c)||c=='('||c==')'||c=='\''||c=='"'||c==';'; }
  void skip(std::size_t depth) {
    for(;;) {
      while(p<s.size() && space(s[p])) ++p;
      if(p<s.size()&&s[p]==';') { while(p<s.size()&&s[p]!='\n') ++p; continue; }
      if(s.substr(p,2)=="#|") {
        p+=2; std::size_t nesting=1;
        while(p<s.size()&&nesting) {
          if(s.substr(p,2)=="#|") { ++nesting; p+=2; }
          else if(s.substr(p,2)=="|#") { --nesting; p+=2; }
          else ++p;
        }
        if(nesting) fail("unterminated block comment");
        continue;
      }
      if(s.substr(p,2)=="#;") { p+=2; read(depth+1); continue; }
      break;
    }
  }
  SExpr read(std::size_t depth) {
    if(depth>limits.max_depth) fail("nesting limit exceeded");
    skip(depth);
    if(p==s.size()) fail("expected expression");
    if(++nodes>limits.max_nodes) fail("node limit exceeded");
    char c=s[p++]; SExpr x;
    if(c=='(') {
      for(;;) {
        skip(depth);
        if(p==s.size()) fail("unterminated list");
        if(s[p]==')') { ++p; return x; }
        if(s[p]=='.' && (p+1==s.size()||delimiter(s[p+1]))) {
          if(x.elements.empty()) fail("dot without list head");
          ++p; x.tail=std::make_shared<SExpr>(read(depth+1)); skip(depth);
          if(p==s.size()||s[p++]!=')') fail("expected close after dotted tail");
          return x;
        }
        x.elements.push_back(read(depth+1));
      }
    }
    if(c==')') fail("unexpected close parenthesis");
    if(c=='\'') {
      if(++nodes>limits.max_nodes) fail("node limit exceeded");
      x.elements.push_back(SExpr::symbol("quote")); x.elements.push_back(read(depth+1)); return x;
    }
    if(c=='"') {
      x.kind=SExpr::Kind::string;
      while(p<s.size()) {
        c=s[p++]; if(c=='"') return x;
        if(c=='\\') {
          if(p==s.size()) fail("unterminated escape");
          c=s[p++];
          switch(c) {
            case 'n': c='\n'; break; case 'r': c='\r'; break; case 't': c='\t'; break;
            case 'a': c='\a'; break; case 'b': c='\b'; break; case 'f': c='\f'; break; case 'v': c='\v'; break;
            case '\n': continue;
            case 'x': {
              unsigned value=0; std::size_t n=0;
              while(p<s.size()&&std::isxdigit(static_cast<unsigned char>(s[p]))) {
                char h=s[p++]; value=value*16+(h<='9'?h-'0':std::tolower(h)-'a'+10);
                if(++n>2) fail("hex escape outside byte range");
              }
              if(n==0||p==s.size()||s[p++]!=';') fail("invalid hex escape");
              c=static_cast<char>(value); break;
            }
            case '\\': case '"': break;
            default: fail("unsupported string escape");
          }
        }
        x.text+=c;
      }
      fail("unterminated string");
    }
    std::size_t start=p-1;
    while(p<s.size()&&!delimiter(s[p])) ++p;
    auto token=s.substr(start,p-start);
    if(token==".") fail("unexpected dot");
    if(token=="#t"||token=="#f") { x.kind=SExpr::Kind::boolean; x.boolean=token=="#t"; return x; }
    auto numeric=token;
    if(numeric.starts_with('+')) numeric.remove_prefix(1);
    std::int64_t i=0;
    auto ir=std::from_chars(numeric.data(),numeric.data()+numeric.size(),i);
    if(ir.ec==std::errc{}&&ir.ptr==numeric.data()+numeric.size()&&!numeric.empty()) { x.kind=SExpr::Kind::integer; x.integer=i; return x; }
    double d=0;
    auto dr=std::from_chars(numeric.data(),numeric.data()+numeric.size(),d);
    if(dr.ec==std::errc{}&&dr.ptr==numeric.data()+numeric.size()&&!numeric.empty()&&std::isfinite(d)) { x.kind=SExpr::Kind::real; x.real=d; return x; }
    const bool number_like=!numeric.empty() && (std::isdigit(static_cast<unsigned char>(numeric[0])) ||
      (numeric[0]=='-'&&numeric.size()>1&&std::isdigit(static_cast<unsigned char>(numeric[1]))));
    if(number_like && (ir.ec==std::errc::result_out_of_range||dr.ec==std::errc::result_out_of_range)) fail("number outside supported range");
    return SExpr::symbol(std::string(token));
  }
public:
  Parser(std::string_view text, ParseLimits lim):s(text),limits(lim) {
    if(s.size()>limits.max_bytes) fail("input size limit exceeded");
  }
  std::vector<SExpr> all() { std::vector<SExpr> result; skip(0); while(p<s.size()) { result.push_back(read(0)); skip(0); } return result; }
};
void write(std::ostream& out,const SExpr& x) {
  switch(x.kind) {
    case SExpr::Kind::symbol: out<<x.text; break;
    case SExpr::Kind::string:
      out<<'"';
      for(unsigned char c:x.text) {
        switch(c) {
          case '\\': out<<"\\\\"; break; case '"': out<<"\\\""; break;
          case '\n': out<<"\\n"; break; case '\r': out<<"\\r"; break; case '\t': out<<"\\t"; break;
          default: if(c<32) out<<"\\x"<<std::hex<<static_cast<unsigned>(c)<<';'<<std::dec; else out<<c;
        }
      }
      out<<'"'; break;
    case SExpr::Kind::integer: out<<x.integer; break;
    case SExpr::Kind::real: {
      std::ostringstream v; v<<std::setprecision(17)<<x.real;
      auto str=v.str(); out<<str; if(str.find_first_of(".eE")==std::string::npos) out<<".0";
      break;
    }
    case SExpr::Kind::boolean: out<<(x.boolean?"#t":"#f"); break;
    case SExpr::Kind::list:
      out<<'('; for(std::size_t i=0;i<x.elements.size();++i) { if(i) out<<' '; write(out,x.elements[i]); }
      if(x.tail) { out<<" . "; write(out,*x.tail); } out<<')'; break;
  }
}
}
std::vector<SExpr> parse_sexprs(std::string_view text, ParseLimits limits) { return Parser(text,limits).all(); }
SExpr parse_sexpr(std::string_view text, ParseLimits limits) {
  auto xs=parse_sexprs(text,limits); if(xs.size()!=1) throw ParseError(0,"expected exactly one expression"); return std::move(xs.front());
}
std::string to_string(const SExpr& x) { std::ostringstream out; write(out,x); return out.str(); }
const SExpr* entry_value(const SExpr& e) {
  if(e.kind!=SExpr::Kind::list) return nullptr;
  if(e.elements.size()==1&&e.tail) return e.tail.get();
  if(e.elements.size()==2&&!e.tail) return &e.elements[1];
  return nullptr;
}
const SExpr* alist_entry(const SExpr& a,std::string_view key) {
  if(a.kind!=SExpr::Kind::list) return nullptr;
  for(const auto& e:a.elements)
    if(e.kind==SExpr::Kind::list&&!e.elements.empty()&&e.elements[0].atom_text()==key && entry_value(e)) return &e;
  return nullptr;
}
const SExpr* alist_value(const SExpr& a,std::string_view key) { auto e=alist_entry(a,key); return e?entry_value(*e):nullptr; }
}
