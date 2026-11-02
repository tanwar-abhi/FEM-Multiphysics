#include "logger.hpp"

Logger::Logger()
{
    // Default constructor
}

Logger::~Logger()
{
    // Default destructor
}

Logger::Logger(const std::string& filePath)
{
    setLogFileName(filePath + logFileName);
}


void Logger::setLogLevel(const std::string& newLogLevel)
{
    std::string newLevel = newLogLevel;
    std::transform(newLevel.begin(), newLevel.end(), newLevel.begin(), ::toupper);
    logLevelStr = newLogLevel;
    if (newLevel == "WARN")
    {
        logLevel = LogLevel::warn;
    }
    else if (newLevel == "DEBUG")
    {
        logLevel = LogLevel::debug;
    }
    else if (newLogLevel == "PROD") {
        logLevel = LogLevel::production;
    }
    else{
        logLevel = LogLevel::info;
        logLevelStr = "INFO";
    }
}

void Logger::setLogFileName(const std::string& fileName)
{
    logFileName = fileName;
}

std::string Logger::getLogFileName() const
{
    return logFileName;
}

std::string Logger::getLogLevelStr() const {
    return logLevelStr;
}


void Logger::printLog(const std::string& logLine)
{
    std::cout<<"["<<logLevelStr<<"] "<<logLine<<std::endl;

}