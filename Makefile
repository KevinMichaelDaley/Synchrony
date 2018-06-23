SRCS := $(wildcard src/*.cpp)
OBJS := $(addprefix tmp_obj/,$(SRCS:%=%.o))
GPU := $(addprefix tmp_gpu/,$(SRCS:%=%.cu.o))
WARNINGS := -Werror -Wall -Wextra -pedantic -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-promo -Wstrict-overflow=5 -Wundef -Wno-unused -Wno-variadic-macros -Wno-parentheses -fdiagnostics-show-option -Wno-error=vla -Wno-unused-command-line-argument -Wno-unused-parameter
all: libsynchrony.so libsynchronyGPU.so

libsynchrony.so: $(OBJS)
	clang-7 -shared $(OBJS) -o $@
libsynchronyGPU.so: $(GPU)
	clang-7 -shared $(GPU) -lcudart_static -L/usr/lib/cuda/lib64 -ldl -lrt -pthread -o $@
tmp_obj/src/%.cpp.o: src/%.cpp
	mkdir -p $(dir $@)
	clang++-7 -c -fPIC -std=c++17 -I./include -g3 $< -o $@ $(WARNINGS)
tmp_gpu/src/%.cpp.cu.o: src/%.cpp
	mkdir -p $(dir $@)
	cp $< tmp_gpu/$<.cu
	clang++-7 -c -fPIC --cuda-gpu-arch=sm_30  -dc -ffast-math -std=c++17 tmp_gpu/$<.cu -o $@ $(WARNINGS) -I./include -O3 --cuda-path=/usr/lib/cuda
clean:
	rm -rf tmp_obj tmp_gpu libsynchrony.so libsynchronyGPU.a
