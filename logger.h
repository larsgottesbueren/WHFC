#pragma once

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>


class Logger {
public:
	explicit Logger(const bool newline) :
			_newline(newline),
			_oss() { }
	template <typename T>
	Logger& operator<< (const T& output) {
		_oss << output << ' ';
		return *this;
	}

	Logger& operator<< (decltype(std::left)& output) {
		_oss << output;
		return *this;
	}

	Logger& operator<< (const decltype(std::setw(1))& output) {
		_oss << output;
		return *this;
	}

	~Logger() {
		std::cout << _oss.str();
		if (_newline) {
			std::cout << std::endl;
		} else {
			std::cout << ' ';
		}
	}

private:
	bool _newline;
	std::ostringstream _oss;
};

class LoggerVoidify {
public:
	void operator& (Logger&) { }
};

#define V(X) #X << "=" << X

#define LOGCC(cond, newline) \
  !(cond) ? (void)0 :        \
  LoggerVoidify() & Logger(newline)

#define LOG  LOGCC(debug, true)
#define LOGWN LOGCC(debug,false)
