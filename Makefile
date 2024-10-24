CC = g++
CFLAGS = -std=c++11
SOURCES = main.cpp linear_eqs.cpp non_linear_eqs.cpp diff_eqs.cpp matrix_inversion.cpp utils.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXEC = numerical_methods_app

all: $(EXEC)

$(EXEC): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(EXEC)
