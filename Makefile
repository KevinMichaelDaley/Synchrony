SRCS := $(wildcard src/*.cpp)
OBJS := $(addprefix tmp_obj/,$(SRCS:%=%.o))
WARNINGS := -Werror -Wall -Wextra -pedantic -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wundef -Wno-unused -Wno-variadic-macros -Wno-parentheses -fdiagnostics-show-option -Wno-error=vla
libsynchrony.so: $(OBJS)
	gcc -shared $(OBJS) -o $@
	
tmp_obj/src/%.cpp.o: src/%.cpp
	mkdir -p $(dir $@)
	g++ -c -fPIC -std=c++17 -I./include -g3 $< -o $@ $(WARNINGS)

clean:
	rm -rf tmp_obj libsynchrony.so
