#include <iostream>
#include "nbsystem.h"

int main(int argc, char *argv[])
{
    NBSystem nbsystem;

    nbsystem.readSettings(argv[1]);

    //std::cout << nbsystem.boundariesString() << std::endl;
    std::cout << nbsystem.stateString() << std::endl;
    return 0;
}
