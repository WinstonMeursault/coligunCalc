# coilgunCalc-L

## 概述

本项目致力于同步感应线圈炮的出口速度与能量转换效率的数值计算

coilgunCalc-L中的"-L"是发射的略写

希望本项目能对尝试数值计算同步感应线圈炮模型的研究者带来少许帮助

本项目的公式依据请参考《电磁感应线圈炮原理与技术》一书^[1]^

## 简要文档

### 开始使用

建议通过

```python
import coilgunCalc-L as cgl
```

来引入coilgunCalc-L

在下面的文档中, 将会以cgl作为coilgunCalc-L的略写

### cgl.basics API

- 线圈自感计算, 定义为:

    ``` python
    def L(coilA, limit=200):
        pass
    ```

    传入线圈对象与limit

    其中limit将传给`scipy.integrate.quad`作为子区间的数量

    增加limit可以使积分结果更加精确,但是达到一定数值后精度没有明显的增加

    建议仅当scipy警告积分结果不够精确时尝试增加这个参数的值

- 线圈互感计算, 定义为:

    ```python
    @lru_cache()
    def M(Ra, Rb, d):
        pass
    ```

    可以通过缓存加速计算(观察矩阵, 由于圆环线圈之间互感只取决于它们的半径与相对位置, 不难得出事实上很多数据是重复的)

    传入两线圈半径(Ra, Rb)与相对位置(d)

- 线圈互感梯度计算, 定义为:

    ```python
    @lru_cache()
    def dM(Ra, Rb, d):
        pass
    ```

    可以通过缓存加速计算(原理与互感计算是相同的)

    传入两线圈半径(Ra, Rb)与相对位置(d)

### cgl.core API

- 驱动线圈类, 定义为:

    ```python
    class drivingCoil():
        def __init__(self, rdi, rde, ld, n, resistivity, Swire, k, s):
            pass

        def R(self):
            #不需要使用者自行调用
            pass
    ```

    `drivingCoil`是驱动线圈的抽象类

    在创建线圈炮抽象类`singleStageCoilgun`与`multiStageCoilgun`均需要传入一个或多个

    下面叙述初始化对象时需要的参数:
        - rdi, rde, ld      驱动线圈的内径/外径/长度
        - n                 驱动线圈匝数
        - resistivity       构成驱动线圈主体的导线材质的电阻率
        - Swire             构成驱动线圈主体的**单根导线**的截面积
        - k                 驱动线圈填充率(定义为截面内导线所占面积与截面总面积之比)

## 参考文献

[1]向红军. 电磁感应线圈炮原理与技术[M].  
