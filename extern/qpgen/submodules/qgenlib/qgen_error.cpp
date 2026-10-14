#include "qgenlib/qgen_error.h"
#include "qgenlib/qgen_except.h"

#include <string>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <inttypes.h>
#include <stdarg.h>

int32_t globalVerbosityThreshold = 100;

void error(const char * msg, ...)
{
  va_list  ap;

  va_start(ap, msg);

  fprintf(stderr, "\nFATAL ERROR - \n");
  vfprintf(stderr, msg, ap);
  fprintf(stderr, "\n\n");

  va_end(ap);

  //throw pexception;
  exit(EXIT_FAILURE);
}

void warning(const char * msg, ...)
{
  va_list  ap;

  va_start(ap, msg);

  time_t current_time;
  char buff[255];
  current_time = time(NULL);

  strftime(buff, 120, "%Y/%m/%d %H:%M:%S", localtime(&current_time));

  fprintf(stderr,"\aWARNING [%s] - ", buff);
  vfprintf(stderr, msg, ap);
  fprintf(stderr,"\n");
  
  va_end(ap);
}

void numerror(const char * msg , ...)
{
  va_list  ap;

  va_start(ap, msg);

  fprintf(stderr,"\nFATAL NUMERIC ERROR - ");
  vfprintf(stderr, msg, ap);
  fprintf(stderr,"\n\n");

  va_end(ap);

  exit(EXIT_FAILURE);
}

void notice(const char * msg, ...) {
  va_list ap;
  va_start(ap, msg);

  time_t current_time;
  char buff[255];
  current_time = time(NULL);

  strftime(buff, 120, "%Y/%m/%d %H:%M:%S", localtime(&current_time));

  fprintf(stderr,"NOTICE [%s] - ", buff);
  vfprintf(stderr, msg, ap);
  fprintf(stderr,"\n");

  va_end(ap);
}

void verbose(int32_t priority, const char * msg, ...) {
  va_list ap;
  va_start(ap, msg);

  if ( globalVerbosityThreshold < priority ) {
    time_t current_time;
    char buff[255];
    current_time = time(NULL);

    strftime(buff, 120, "%Y/%m/%d %H:%M:%S", localtime(&current_time));

    fprintf(stderr,"VERBOSE_MESSAGE_%d [%s] - ", priority, buff);
    vfprintf(stderr, msg, ap);
    fprintf(stderr,"\n");
  }

  va_end(ap);
}

void catprintf(std::string &s, const char * msg, ...)
{
  va_list ap;

  va_start(ap, msg);

  char buf[1000];
  vsnprintf(buf, 1000, msg, ap);

  s += buf;
  va_end(ap);
}

int32_t cat_join_int32(std::string& s, std::vector<int32_t>& v, const char* delim) {
  char buf[65535];
  int32_t len = 0, offset = 0;
  int32_t n = (int32_t)v.size();
  if ( n > 0 ) {
    offset = snprintf(buf, 65535, "%d", v[0]);
    for(int32_t i=1; i < n; ++i) {
      offset += snprintf(buf + offset, 65535, "%s%d", delim, v[i]);
      if ( offset > 65000 ) { // copy the string and reset
        len += offset;
        s += buf;
        buf[0] = '\0';
        offset = 0;
      }
    }
    s += buf;
  }
  return len + offset;
}

int32_t cat_join_uint64(std::string& s, std::vector<uint64_t>& v, const char* delim) {
  char buf[65535];
  int32_t len = 0, offset = 0;
  int32_t n = (int32_t)v.size();
  if ( n > 0 ) {
    offset = snprintf(buf, 65535, "%" PRIu64, v[0]);
    for(int32_t i=1; i < n; ++i) {
      offset += snprintf(buf + offset, 65535, "%s%" PRIu64, delim, v[i]);
      if ( offset > 65000 ) { // copy the string and reset
        len += offset;        
        s += buf;
        buf[0] = '\0';
        offset = 0;
      }
    }
    s += buf;
  }
  return len + offset;  
}

int32_t cat_join_str(std::string& s, std::vector<std::string>& v, const char* delim) {
  char buf[65535];
  int32_t len = 0, offset = 0;
  int32_t n = (int32_t)v.size();
  if ( n > 0 ) {
    offset = snprintf(buf, 65535, "%s", v[0].c_str());
    for(int32_t i=1; i < n; ++i) {
      if ( offset + v[i].size() > 65536 ) { // copy the string and reset
        len += offset;
        s += buf;
        buf[0] = '\0';
        offset = 0;
      }      
      offset += snprintf(buf + offset, 65535, "%s%s", delim, v[i].c_str());
    }
    s += buf;
  }
  return offset + len;
}
