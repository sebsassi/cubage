CXX = g++
CXXFLAGS = -Wall -Wextra -Wpedantic -Wconversion -Werror -std=c++20 -O3 -pg -g

BENCH = benchmark_genz_malik

BENCH_SUFFIX = .bench

$(BENCH): %: %.cpp
	$(CXX) $(CXXFLAGS) -o $@$(BENCH_SUFFIX) $<

benchmarks: $(BENCH)

clean:
	rm $(BENCH)