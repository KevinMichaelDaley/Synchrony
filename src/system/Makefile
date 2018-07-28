WARNINGS := -Werror -Wall -Wextra -pedantic -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-promo -Wstrict-overflow=5 -Wundef -Wno-unused -Wno-variadic-macros -Wno-parentheses -fdiagnostics-show-option -Wno-error=vla -Wno-error=vla-extension -Wno-unused-command-line-argument -Wno-unused-parameter
GPU_SRCS := $(wildcard ../*.cpp)
gpu:
	clang++-6.0 -Xcuda-ptxas "-v" -DNDEBUG -O3 --cuda-gpu-arch=sm_61  -ffast-math -std=c++17 $(WARNINGS) -I ../../include -I/usr/local/cuda-8.0/include  TacomaNarrows.cu -x cuda  --cuda-path=/usr/local/cuda-8.0 -o tacoma -g3 -O3 -L/usr/local/cuda-8.0/lib64 -lcudart_static -ldl -lrt -pthread 
