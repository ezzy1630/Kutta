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
#include <sys/stat.h>
#include <unistd.h>
#endif

namespace kutta::app {
// No worker thread and no blocking stdin read. Handlers only set a flag; the
// solver observes it between outer iterations and writes normal final outputs.
// EOF means "stop" only on a pipe, which is how the app holds the sidecar's
// stdin. A file, /dev/null, a terminal or no stdin at all would otherwise turn
// a scripted run into an immediate zero-iteration stop.
class StopRequest {
 public:
  explicit StopRequest(bool watch_stdin) : watch_stdin_(watch_stdin && stdin_is_pipe()) {
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
    DWORD available=0;
    if(!PeekNamedPipe(input,nullptr,0,nullptr,&available,nullptr))return true;
    if(available){char bytes[4096];DWORD read=0;if(!ReadFile(input,bytes,available<sizeof(bytes)?available:DWORD(sizeof(bytes)),&read,nullptr))return true;}
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
  static bool stdin_is_pipe() {
#ifdef _WIN32
    const auto input=GetStdHandle(STD_INPUT_HANDLE);
    return input!=INVALID_HANDLE_VALUE && input!=nullptr && GetFileType(input)==FILE_TYPE_PIPE;
#else
    struct stat s{};
    return fstat(STDIN_FILENO,&s)==0 && (S_ISFIFO(s.st_mode) || S_ISSOCK(s.st_mode));
#endif
  }
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
