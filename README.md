 # In progress...

# Sun-Earth/Moon L1 Lissajous orbit multiple shooting

## CRTBP model

Third order approximation around L1 point --> CRTBP synodic frame (done)

![](Test_02_CRTBP.gif)

（分了400段，各点位置精度1e-9，总速度精度1e-7，然后level 2用17次打靶完成，总共128秒）

## Ephemeris model


![](Test_02_Ephemeris.gif)

（日地月星历下，同样的分了400段，各点位置精度10m，总速度精度0.1m/s，然后level 2用4次打靶完成，总共200秒，很快啊，不过这次开了4核8线程并行。）




# Plan 
- ode45 和 ode113 的算法比较来验证计算精度的设计
- add basic solar radiation model, maybe
- automatically parse some critical input parameters
- more testing examples
- improve performance by benchmarking the performance
