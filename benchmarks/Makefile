CXX = g++
CXXFLAGS = -Wall -Wextra -Wpedantic -Wconversion -Werror -std=c++20 -O3 -pg -g

BENCH = benchmark_genz_malik

$(BENCH): %: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

benchmarks: $(BENCH)

clean:
	rm $(BENCH)