#include <fstream>
#include <iostream>
#include <string>
#pragma once


class FASTAParse {
public:
    std::string mHeader;
    std::string mSequence;
    FASTAParse(const std::string &inputDir);
};
