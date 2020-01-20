 
#ifndef _COMMAND_H
#define _COMMAND_H

#include <math.h>

class Command {
    public:
        virtual ~Command() {};
        virtual int main(int argc, char *argv[]) = 0;
    };

#endif // _COMMAND_H
