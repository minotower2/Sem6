QT += widgets
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3
QMAKE_CXXFLAGS += -mfpmath=sse -fstack-protector-all -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format -std=c++17

HEADERS = window.h functions.h newton_approximation.h spline_approximation.h
SOURCES = main.cpp window.cpp functions.cpp newton_approximation.cpp spline_approximation.cpp

zip.commands = zip graph.zip $$SOURCES $$HEADERS graph.pro
QMAKE_EXTRA_TARGETS += zip