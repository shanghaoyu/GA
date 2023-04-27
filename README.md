# GA
这个程序使用遗传进化算法进行优化核力参数 
# gasrc
其中gasrc是遗传进化算法的源文件 
* 使用了openmp并行，通过make进行编译出“ga”，在可并行的节点上注意调整一下静态omp-set
* “ga”需要两个个输入文件“constan.d”“expdata.d"和一个可执行文件“phase”
* “constan.d”中给出了优化参数空间的维数，和参数的上下界，注意这个维数需要和所调用“phase”的参数个数一致
* “expdata.d"前面部分指出拟合哪些相移
* “phase”计算相移，注意其顺序应与“expdata.d"一致
# phasesrc
phasesrc给出了一个可以使用phase可执行文件的实例
