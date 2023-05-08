CXX=g++
CXXFLAGS=-std=c++17 -flto -O3 -march=native

LD=g++
LDFLAGS=-flto=auto

SRCSDIR=src
OBJSDIR=obj
DESTDIR=target

SRCS=$(wildcard $(SRCSDIR)/*.cpp)
OBJS=$(patsubst $(SRCSDIR)/%.cpp,$(OBJSDIR)/%.o,$(SRCS))
DEPS=$(patsubst $(SRCSDIR)/%.cpp,$(OBJSDIR)/%.d,$(SRCS))

RIblast: $(OBJS)
	@mkdir -p $(DESTDIR)
	$(LD) $(LDFLAGS) -o $(DESTDIR)/$@ $^ $(LDLIBS)

-include $(DEPS)

$(OBJSDIR)/%.o: $(SRCSDIR)/%.cpp
	@mkdir -p $(OBJSDIR)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

c: clean
clean:
	rm -rf $(OBJSDIR) $(DESTDIR)
