# Command Line Arguments using qgenlib

## Introduction

`qgenlib` provides a macro-based command line parsing library. 
This library was first implemented by [Gonçalo Abecasis](https://en.wikipedia.org/wiki/Gon%C3%A7alo_Abecasis) and used in many [tools developed by Center for Statistical Genetics](https://github.com/statgen/libStatGen) at the 
University of Michigan. 

The key features of the command line parsing library is that the user can define the command line options using a simple macro and the library will automatically generate the help message and parse the command line arguments.

## Example Usage

The following is a simple example to demonstrate how to use the command line parsing library.

```cpp
#include <vector>
#include "qgenlib/params.h"
int main(int argc, char** argv) {
    // define the variables to store the command line arguments
    std::string str1;
    std::vector<std::string> str2;
    int int3 = 0;
    double double4 = 0.0;
    bool bool5 = false;

    // the macro definition starts with defining a paramList object.
    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
        LONG_PARAM_GROUP("Options for string arguments", NULL)
        LONG_STRING_PARAM("str1", &str1, "A dummy string argument")
        LONG_MULTI_STRING_PARAM("str2", &str2, "Another dummy string argument (can be specified multiple times)")

        LONG_PARAM_GROUP("Options for non-string arguments", NULL)
        LONG_INT_PARAM("int3", &int3, "A dummy integer argument")    
        LONG_DOUBLE_PARAM("double4", &double4, "A dummy double argument")
        LONG_PARAM("bool5",&bool5,"A dummy boolean argument")
    END_LONG_PARAMS();

    // 
    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    printf("str1: %s\n", str1.c_str());
    printf("str2: ");
    for (int i = 0; i < str2.size(); i++) {
        printf("%s ", str2[i].c_str());
    }
    printf("\n");
    printf("int3: %d\n", int3);
    printf("double4: %f\n", double4);
    printf("bool5: %s\n", bool5 ? "true" : "false");

    return 0;
}
```

This example can also be found in the [qbcf](https://github.com/hyunminkang/qbcf/) repository in the `qbcf test-qgenlib` tool. 

## Expected Output 

Output of an example run with some parameter settings are shown as follows:

```sh
qbcf test-qgenlib --str1 10 --str2 abc --str2 def --int3 12345 --bool5
[qbcf test-qgenlib] -- Test basic functionalities in qgenlib

 Copyright (c) 2009-2022 by Hyun Min Kang and Adrian Tan
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Available Options

The following parameters are available. Ones with "[]" are in effect:
       Options for string arguments : --str1 [10], --str2 [abc, def]
   Options for non-string arguments : --int3 [12345], --double4 [0.00],
                                      --bool5 [ON]

Run with --help for more detailed help messages of each argument.

str1: 10
str2: abc def 
int3: 12345
double4: 0.000000
bool5: true
```

Printing out the help message can be done with `--help` option as follows:

```sh
qbcf test-qgenlib --help                                              
[qbcf test-qgenlib] -- Test basic functionalities in qgenlib

 Copyright (c) 2009-2022 by Hyun Min Kang and Adrian Tan
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Options for string arguments ==
   --str1    [STR: ]             : A dummy string argument
   --str2    [V_STR: ]           : Another dummy string argument (can be specified multiple times)

== Options for non-string arguments ==
   --int3    [INT: 0]            : A dummy integer argument
   --double4 [FLT: 0.00]         : A dummy double argument
   --bool5   [FLG: OFF]          : A dummy boolean argument


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```

