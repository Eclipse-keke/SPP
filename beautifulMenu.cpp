#include"myhead.h"

void printCentered(const std::string& text, int width)
{
    int padding = (width - text.size()) / 2;
    for (int i = 0; i < padding; ++i)
        std::cout << ' ';
    std::cout << text << std::endl;
}

void printBorderLine(int length)
{
    std::cout << "¨X";
    for (int i = 0; i < length - 2; ++i)
        std::cout << "¨T";
    std::cout << "¨[" << std::endl;
}

void printMiddleBorderLine(int length)
{
    std::cout << "¨c";
    for (int i = 0; i < length - 2; ++i)
        std::cout << "©¤";
    std::cout << "¨f" << std::endl;
}

void printBottomBorderLine(int length)
{
    std::cout << "¨^";
    for (int i = 0; i < length - 2; ++i)
        std::cout << "¨T";
    std::cout << "¨a" << std::endl;
}

void printMenuItem(const std::string& item, int width)
{
    int padding = (width - item.size() - 2) / 2;
    std::cout << "¨U";
    for (int i = 0; i < padding; ++i)
        std::cout << ' ';
    std::cout << item;
    for (int i = 0; i < width - item.size() - 2 - padding; ++i)
        std::cout << ' ';
    std::cout << "¨U" << std::endl;
}