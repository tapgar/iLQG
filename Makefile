CC	:= g++
LD	:= $(CC)
MJDIR	:= /home/tapgar/libcassie-master/mjpro150

INC	:= -Iinclude -I$(MJDIR)/include -I/home/tapgar/eigen
CFLAGS	:= -std=c++11 -Wall
LDFLAGS	:= -L$(MJDIR)/bin -Wl,-rpath,'$(MJDIR)/bin'
LIBS	:= -lmujoco150 -lglew -lGL -l:libglfw.so.3

SRC	:= $(wildcard src/*.cpp)
OBJ	:= $(patsubst src/%.cpp,obj/%.o,$(SRC))

vpath %.cpp src

define make-goal
$1/%.o: %.cpp
	$(CC) $(CFLAGS) $(INC) -MMD -c $$< -o $$@
endef

all: checkdirs build

clean:
	rm -rf obj/
	rm -f test

DIP: 
	$(CC) $(CFLAGS) main_dip.cpp $(SRC) $(INC) $(LDFLAGS) $(LIBS) -DnQ=3 -DnQd=3 -DnU=1 -o DIPilqg

Hopper: 
	$(CC) $(CFLAGS) main_hopper.cpp $(SRC) $(INC) $(LDFLAGS) $(LIBS) -DnQ=6 -DnQd=6 -DnU=3 -DnL=3 -DnP=150 -o Hopperilqg

Cassie: 
	$(CC) $(CFLAGS) main_cassie.cpp $(SRC) $(INC) $(LDFLAGS) $(LIBS) -DnQ=32 -DnQd=32 -DnU=10 -DnL=10 -DnP=500 -o Cassieilqg

SLIP: 
	$(CC) $(CFLAGS) main_slip.cpp $(SRC) $(INC) $(LDFLAGS) $(LIBS) -DnQ=9 -DnQd=9 -DnU=4 -DnL=10 -DnP=500 -o SLIPilqg

build: $(OBJ)
	$(LD) $^ -o $(OUT) $(LDFLAGS) $(LIBS)

checkdirs: obj

obj:
	@mkdir -p $@

$(foreach bdir,obj,$(eval $(call make-goal,$(bdir))))

.PHONY: all checkdirs clean

-include $(OBJ:%.o=%.d)
