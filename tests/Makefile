CXX = g++
CXXFLAGS = -Wall -Wextra -Wpedantic -Wconversion -Werror -std=c++20 -Og -g

TESTS = test_cubage test_gauss_kronrod test_genz_malik

$(TESTS): %: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

tests: $(TESTS)

clean:
	rm $(TESTS)