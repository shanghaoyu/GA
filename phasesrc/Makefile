# 设置编译器
FC = ifort

# 设置源文件
SRCS = add.f main.f n2loadd.f phase.f

# 设置输出文件名
OUT = phase

# 编译规则
all: $(OUT)

$(OUT): $(SRCS)
	$(FC) $(SRCS) -o $(OUT)

clean:
	rm -f $(OUT) *.mod