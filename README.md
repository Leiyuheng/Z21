# Z21  
主要参考文献：  
杨儒贵——《高等电磁理论》  
Chen Yongpin —— 《A New Green’s Function Formulation for Modeling Homogeneous Objects in Layered Medium》  
ALIREZA H.MOHAMMADIAN —— 《A Theoretical and Experimental Study of Mutual Coupling in Microstrip Antenna Arrays》  
K. Sarabandi —— 并矢分析相关内容  
  
基于空腔模型与反应原理计算两个贴片天线之间的互耦合  
选用的格林函数为半自由空间格林函数  
贴片天线的工作模式设定为TM01模式  
%%====================首次上传记录==========================  
  
目前计算的磁流强度为sin(kx)  
并非严格的空腔模型解  
  
文件具体说明： 
3D文件中是计算电流源产生的电场与计算磁流源产生的磁场（理想点源） 
分别对应E和H两个文件夹  
文件中主要对应三个部分 分别是主函数、格林函数计算部分与电场/磁场计算部分  
文件夹中main_single是指计算单个理想点源的电场/磁场  
main_rotate是指计算线电流在不同角度下的电场 使用的方法是求离散点后插值进行计算（但这种方法求解出的数量级不对，因为剖分问题，不建议精确计算使用，可以查看变化趋势，总体来说不是很准）  

而half_space_gf与half_space_gf_cal是使用两种不同方式计算的格林函数，起计算结果应当是完全一样的  
唯一的区别在于half_space_gf是直接使用矩阵进行计算，如果要计算积分则涉及到剖分问题  
所以后续修改了代码 使用syms进行编写 便有了half_space_gf_cal这一版本代码  

E_calculate与H_calculate是计算点源电场/磁场的函数 较为简单

而文件中的main与mutual_calculate是结合了上述修改后的代码，计算两条线电流之间的互耦合（后言：这一版本中mutual_calculate关于GF中yz与xz分量处理有误,46和49行应该是负号）  

%% ==========================第二次上传记录=============================
  
更新后的的GF是计算理想磁流之间的互阻抗  
并添加了查看磁流旋转的部分  
并修正了3D文件中的部分错误  
两张图片便是线磁流平行或共线是旋转其中一条获得的阻抗曲线  