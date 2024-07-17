#ifndef TIMER_HPP
#define TIMER_HPP

#include <iostream>
#include <chrono>
#include <string>



class Timer
{
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    std::chrono::duration<double> duration;
public:

    std::string timerName;

    Timer(std::string timerName1 = "Timer");
    ~Timer();

    void stopPrint();
};

// create timer and set starting time. Optionally set Timer name (default "Timer")
Timer::Timer(std::string timer_name) {
        timerName = timer_name;
        start = std::chrono::high_resolution_clock::now();
}

Timer::~Timer() {
    stopPrint();
}

// stop and print the time passed since object created (0.0002ms lost in creation)
void Timer::stopPrint() {
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;

    std::cout << timerName << " - ";
    
    std::cout << duration.count() * 1000 << "ms \n";
}

#endif