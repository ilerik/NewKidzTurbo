#ifndef TURBO_Utility_Timer
#define TURBO_Utility_Timer

#include <cassert>
#include <chrono>
#include <mutex>
#include <memory>
#include <thread>

//Simple timer object
class Timer {

	std::chrono::time_point<std::chrono::high_resolution_clock> _pauseTime;
public:
	bool IsActive;
	std::chrono::time_point<std::chrono::high_resolution_clock> StartTime;
	std::chrono::time_point<std::chrono::high_resolution_clock> EndTime;
	std::chrono::high_resolution_clock::duration ElapsedTime;

	Timer() {
		IsActive = false;
		ElapsedTime = std::chrono::high_resolution_clock::duration(0);
	};

	double ElapsedTimeMilliseconds() {
		return std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(ElapsedTime).count();
	};

	void Start() {
		StartTime = std::chrono::system_clock::now();
		_pauseTime = StartTime;
		ElapsedTime = std::chrono::high_resolution_clock::duration(0);
	};

	void Resume() {
		if (!IsActive) {
			_pauseTime = std::chrono::system_clock::now();
			IsActive = true;
		};
	};

	void Pause() {
		if (IsActive) {
			ElapsedTime += std::chrono::system_clock::now() - _pauseTime;
			IsActive = false;
		};
	};

	void Reset() {
		if (!IsActive) {
			ElapsedTime = std::chrono::high_resolution_clock::duration(0);
		};
	};

	void Stop() {
		Pause();
		EndTime = std::chrono::system_clock::now();
	};
};

#endif