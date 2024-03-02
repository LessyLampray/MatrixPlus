GCC = g++ -std=c++17 -Wall -Wextra -Werror 
LIBS := -lgtest
OS := $(shell uname -s)
CFLAGS := -ggdb -lstdc++ -pthread

ifeq ($(OS), Darwin)
 PROJECT_NAME = s21_matrix_oop.a
endif

ifeq ($(OS), Linux)
 PROJECT_NAME = s21_matrix_oop.cc
 LIBS += $(CFLAGS)
endif

all: test

test: clean s21_matrix_oop.a
	@echo "┏=========================================┓"
	@echo "┃     Testing ...                         ┃"
	@echo "┗=========================================┛"
	@$(GCC) $(PROJECT_NAME) test.cc -o test $(LIBS) 
	@./test

valgrind: test
	@valgrind --tool=memcheck --leak-check=yes  ./test

s21_matrix_oop.a:
	@$(GCC) -c s21_*.cc
	@ar rc s21_matrix_oop.a s21_*.o
	@ranlib s21_matrix_oop.a
	@rm *.o

leak: clean
	@echo "┏=========================================┓"
	@echo "┃     Checking for memory leaks ...       ┃"
	@echo "┗=========================================┛"
	@$(GCC) test.cc s21_*.cc -lgtest -o leaks -g
	@CK_FORK=no leaks --atExit -- ./leaks
	@rm ./leaks

style_check:
	@echo "┏=========================================┓"
	@echo "┃  Checking your code for Google Style    ┃"
	@echo "┗=========================================┛"
	@clang-format -style=Google -n *.cc *.h

format_style:
	@echo "┏=========================================┓"
	@echo "┃  Formatting your code for Google Style  ┃"
	@echo "┗=========================================┛"
	@find . \( -name '*.h' -o -name '*.cc' \) -print0 | xargs -0 clang-format -i -style=Google

clean:
	@rm -rf s21_*o
	@rm -rf *.a
	@rm -rf a.out
	@rm -rf test.*o
	@rm -rf report
	@rm -rf test
	@rm -rf *.gcda
	@rm -rf *.dSYM
	@rm -rf leaks

.PHONY: all test valgrind s21_matrix_oop.a gcov_report leak clang clean