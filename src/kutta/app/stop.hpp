#pragma once

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <cerrno>
#include <csignal>
#include <poll.h>
#include <unistd.h>
#endif

namespace kutta::app {
// No worker thread and no blocking stdin read. Handlers only set a flag; the
// solver observes it between outer iterations and writes normal final outputs.
class StopRequest {
 public:
  explicit StopRequest(bool watch_stdin) : watch_stdin_(watch_stdin) {
#ifdef _WIN32
    InterlockedExchange(&requested_,0);
    installed_=SetConsoleCtrlHandler(handler,TRUE)!=0;
#else
    requested_=0;
    struct sigaction action{};
    action.sa_handler=handler;sigemptyset(&action.sa_mask);
    installed_=sigaction(SIGINT,&action,&previous_)==0;
#endif
  }
  ~StopRequest() {
#ifdef _WIN32
    if(installed_) SetConsoleCtrlHandler(handler,FALSE);
#else
    if(installed_) sigaction(SIGINT,&previous_,nullptr);
#endif
  }
  StopRequest(const StopRequest&)=delete;
  StopRequest& operator=(const StopRequest&)=delete;
  bool requested() {
#ifdef _WIN32
    if(InterlockedCompareExchange(&requested_,0,0))return true;
    if(!watch_stdin_)return false;
    const auto input=GetStdHandle(STD_INPUT_HANDLE);
    const auto type=GetFileType(input);
    if(type==FILE_TYPE_PIPE) {
      DWORD available=0;
      if(!PeekNamedPipe(input,nullptr,0,nullptr,&available,nullptr))return true;
      if(available){char bytes[4096];DWORD read=0;if(!ReadFile(input,bytes,available<sizeof(bytes)?available:DWORD(sizeof(bytes)),&read,nullptr))return true;}
    } else if(type==FILE_TYPE_DISK) {
      LARGE_INTEGER zero{},position{},size{};
      if(GetFileSizeEx(input,&size) && SetFilePointerEx(input,zero,&position,FILE_CURRENT))return position.QuadPart>=size.QuadPart;
    } else if(type==FILE_TYPE_CHAR) {
      DWORD mode=0;
      if(!GetConsoleMode(input,&mode))return true; // NUL, rather than a console
    } else if(type==FILE_TYPE_UNKNOWN || input==INVALID_HANDLE_VALUE || input==nullptr)return true;
    return false;
#else
    if(requested_)return true;
    if(!watch_stdin_)return false;
    pollfd input{STDIN_FILENO,POLLIN|POLLHUP,0};
    const auto ready=poll(&input,1,0);
    if(ready<0)return errno!=EINTR;
    if(!ready)return false;
    if(input.revents&(POLLHUP|POLLERR|POLLNVAL))return true;
    if(input.revents&POLLIN){char bytes[4096];const auto n=read(STDIN_FILENO,bytes,sizeof(bytes));return n==0 || (n<0 && errno!=EINTR && errno!=EAGAIN);}
    return false;
#endif
  }
 private:
  bool watch_stdin_,installed_=false;
#ifdef _WIN32
  inline static volatile LONG requested_=0;
  static BOOL WINAPI handler(DWORD event) {
    if(event!=CTRL_C_EVENT && event!=CTRL_BREAK_EVENT)return FALSE;
    InterlockedExchange(&requested_,1);return TRUE;
  }
#else
  inline static volatile std::sig_atomic_t requested_=0;
  struct sigaction previous_{};
  static void handler(int) {requested_=1;}
#endif
};
} // namespace kutta::app
