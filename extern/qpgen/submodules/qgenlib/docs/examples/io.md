# I/O Functions

## Introduction

`qgenlib` provides functions for basic input/output functions to handle messages to print out during the program execution.

## Error, Warning, and Notification Messages

Users can output `printf`-style functions, with variable number of arguments, to output error, warning, and notification messages. This function acts like `fprintf` in the sense that the message is printed to the standard error output, with the following differences:

* `notice()` function automatically prints the timestamp before each message.
* `warning()` function automatically prints the timestamp, give a beep sound, and prints and the word "WARNING" before each message.
* `error()` function prints out the error message and terminate the program.

Below are examples of how to use these functions:

```cpp
// int n : number of lines processed
notice("Processed %d lines of input lines\n", n);
// Example output: 
// NOTICE [2024/02/29 14:21:09] - Processed 1000 lines of input lines

// int n : current line number
warning("Line %d is not properly formatted\n", n);
// Example output: (with beep)
// WARNING [2024/02/29 14:21:09] - Line 1234 is not properly formatted

// int n : current line number
error("Cannot parse line %d\n", n);
// Example output: (exits the program)
// FATAL ERROR - 
// Cannot parse line 1234
```

Some of these examples can also be found in the [qbcf](https://github.com/hyunminkang/qbcf/) repository in the `qbcf test-qgenlib` tool.

## Other `printf`-style Functions

`qgenlib` also provides other printf-style functions to handle messages. 

* `catprintf()` function appends printf-style string to the current string.
* `hprintf()` function write messages to `htsFile*` object.

```cpp
// append a string to a std::string object in print-style
std::string str("Hello");
int num = 1234;
std::string name("Customer");
// str will contains "Hello, Customer No. 1234" after running the following line
catprintf(str, ", %s No. %d", name.c_str(), num); 

// write a message to a htsFile* object in printf-style
htsFile* fp = hts_open("output.txt", "w");
hprintf(fp, "Hello, %s No. %d", name.c_str(), num); 
hts_close(fp);
```