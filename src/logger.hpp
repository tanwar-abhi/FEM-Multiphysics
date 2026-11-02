#pragma once

#include <iostream>
#include <string>
#include <algorithm>

enum class LogLevel {
    info,
    warn,
    debug,
    production
};



class Logger 
{
    private:
    LogLevel logLevel = LogLevel::info;
    std::string logLevelStr = "INFO";

    std::string logFileName = "logFile.txt";

    public:
    Logger();
    Logger(std::string& filePath);
    ~Logger();

    Logger(const std::string &level);

    void setLogLevel(const std::string& level);

    // LogLevel getLogLevel();
    std::string getLogLevelStr() const;

    void setLogFileName(const std::string& name);
    std::string getLogFileName() const;

    void printLog(const std::string& str);
};