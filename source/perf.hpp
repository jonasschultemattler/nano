#include <sys/resource.h>

class PerformanceTracker {
private:
    std::chrono::steady_clock::time_point startTime;
    size_t initialMemory;

public:
    void start() {
        startTime = std::chrono::steady_clock::now();
        initialMemory = getCurrentMemoryUsage();
    }

    void report() {
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
            (std::chrono::steady_clock::now() - startTime);

        size_t currentMemory = getCurrentMemoryUsage();

        std::cout << "Execution Time: " << duration.count() << "ms" << std::endl;
        std::cout << "Memory Used: " << (currentMemory - initialMemory) << " bytes" << std::endl;
    }

    size_t getCurrentMemoryUsage() {
        // Platform-specific memory retrieval
        // Implementation varies by system
        struct rusage usage;
        int ret;
        ret = getrusage(RUSAGE_SELF, &usage);
        return usage.ru_maxrss;
    }
};