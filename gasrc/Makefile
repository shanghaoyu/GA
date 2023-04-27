# 设置编译器
FC = ifort

FFLAGS = -qopenmp

# 设置其他源文件
SRCS =constant_module.f phase_module.f fun.f ga.f pikaia.f 

# 设置输出文件名
OUT = ga

# 编译规则
all: $(OUT)

$(OUT): $(SRCS)
	$(FC) $(FFLAGS) $(SRCS) -o $(OUT)

clean:
	rm -f $(OUT)
