<script type="text/javascript" async src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-MML-AM_CHTML"> </script>
formula1: $$n==x$$

formula2: $$n!=x$$

formula3: (m==y)

formula4: [m!=y]

formula5: \(k==z\)

formula6: \[k!=z\]

<table><tr><td bgcolor=orange><font face= "黑体" > 一、迭代法</td></tr></table>

引言：在学习CFD中非常关键的一步是了解离散后的代数方程的求解方法，这也是连接学习CFD理论和编写CFD程序之间的纽带. 如果能在习得迭代法后再学习CFD理论，可在学习理论的同时完成程序的编写，非常有成就感. 本文将介绍最简单的几种迭代法，对于编写小程序练手已经够用.

对于中小型方程组，一般可采用直接法求解，如果不考虑舍入误差，直接法能在预定的计算步骤内得到方程的精确解，计算效率高且可靠.而对于维数较高，尤其是大型的稀疏矩阵（0元素较多）方程组，在CFD中碰到矩阵多为稀疏矩阵，由于直接法计算量太大，因此可采用迭代法进行.常见的迭代法有<font face="times new roman">Jacobi</font>迭代法，<font face="times new roman">Gauss-Seidel</font>迭代法，基本迭代法的加速方法（<font face="times new roman">SOR</font>）共轭梯度法，最速下降法等，对于非线性方程（组）的迭代法，本文将不涉及.

准备：对于线性方程组：
$$
Ax=b\tag{1}
$$

其中$A\in R^ {n\times n}$, $f\in R^{n\times n}$, $x\in R^{n}$ 

我们可以将其转化为等价方程：

$$
x=Bx+f\tag{2}
$$

进一步的：
$$
x^{(k+1)}=Bx^{(k)}+f
\tag{3}
$$


其中，$B\in R^{n\times n}$, $f\in R^n$, $x\in R^n$,矩阵$B$即为迭代矩阵，不同的迭代方法的区别在于迭代矩阵$B$的不同.可以看出求解方程组解的过程就是求解向量序列极限的过程，极限存在则迭代法收敛，否则迭代发散.

其一般迭代过程为：

$$
x^{(k+1)}=Bx^{(k)}+f\tag{4}
$$

---

<font face= "times new roman" color=blue > 1.1 Jacobi迭代法 </font>
$$
a_{11}x_1+a_{12}x_2+\cdots+a_{1n}x_n=b_1\\

a_{21}x_1+a_{22}x_2+\cdots+a_{2n}x_n=b_2\\

\cdots\\

a_{n1}x_1+a_{n2}x_2+\cdots+a_{nn}x_n=b_n\\
\tag{5}
$$

其中矩阵$A=(a_{ij})_{n\times n}$为非奇异矩阵，且$a_{ii}\neq 0$

变形后：

$$
x_1=\frac{1}{a_{11}}(b_1-\sum_{j=2}^na_{1j}x_j)\\

x_2=\frac{1}{a_{22}}(b_2-\sum_{j=1,j\neq2}^na_{1j}x_j)\\

\cdots\\

x_i=\frac{1}{a_{ii}}(b_i-\sum_{j=1,j\neq i}^na_{ij}x_j)\\

\cdots\\

x_n=\frac{1}{a_{nn}}(b_n-\sum_{j=1,j\neq n}^{n}a_{nj}x_j)
\tag{6}
$$

可以将其写出迭代格式：

$$
x_i^{(k+1)}=x^k+\frac{1}{a_{ii}}(b_i-\sum_{j=1}^na_{ij}x_j^{(k)}),(i=1,2,\cdots,n)
\tag{7}
$$
进一步的可得到其迭代矩阵：


$$
B_J=D^{-1}(L+U),f=D^{-1}b\tag{8}
$$

其中$A=D-L-U，D=diag(a_{11},a_{22},\cdots,a_{nn}),L$为$A$矩阵下三角元素的负数，$U$为$A​$矩阵上三角元素的负数.

例：用<font face="times new roman">Jacobi迭代法求解下述方程组</font>
$$
\begin{bmatrix}
8&-3&2\\4&11&-1\\2&1&4
\end{bmatrix}\begin{bmatrix}x_1\\x_2\\x_3\end{bmatrix}=\begin{bmatrix}20\\33\\12\end{bmatrix}\tag{9}
$$

该方程有精确解$X=(3,2,1)$

使用fortran程序求解，代码如下：

```fortran
module jacobi
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010-4-5
!-----------------------------------------------------
!  Description : jacobi迭代法模块
!    
!-----------------------------------------------------
!  Parameters  :
!      1.    IMAX最大允许迭代次数  
!      2.    tol误差容限
!  Contains    :
!      1.    solve 雅克比迭代法方法函数
!      2.
!-----------------------------------------------------
!  Post Script :
!      1.
!      2. 
!-----------------------------------------------------

implicit real*8(a-z)
integer::IMAX=200
real*8::tol=1d-7

contains

subroutine solve(A,b,x,x0,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  雅克比迭代法函数
!               用于计算方程 AX=b
!-----------------------------------------------------
!  Input  parameters  :
!       1.  A,b 意义即  AX=b
!       2.  x0迭代初值
!       3.  N 方程的维数
!  Output parameters  :
!       1. x 方程的解
!       2.
!  Common parameters  :
!
!----------------------------------------------------

implicit real*8(a-z)
integer::N
integer::i,j,k

real*8::A(N,N),b(N),x(N),x0(N)

real*8::x1(N),x2(N)

x1=x0

!写入标题
  write(102,501)
  501 format(//,18x,'Jacobi迭代法',//)

do k=1,IMAX

   do i=1,N
        s=0
        
        do j=1,N
          if (j==i) cycle           
          s=s+A(i,j)*x1(j)      
        
        end do
       
       x2(i)=(b(i)-s)/A(i,i)  

   end do


! 这段程序用于判断精度，满足精度时退出循环，使用向量的二范数进行精度的判断.
   dx2=0
   do i=1,N
    dx2=dx2+(x1(i)-x2(i))**2
   end do
   dx2=dsqrt(dx2)

   if (dx2<tol)  exit
!----------------------------------------      
   x1=x2

  !记录迭代中间值
     write(102,502)k,x1
     502 format(I3,3F12.8)
  !----
end do

x=x2

end subroutine solve


end module jacobi



program  main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-5
!-----------------------------------------------------
!  Purpose   :  采用雅克比迭代法计算线性方程
!    
!-----------------------------------------------------
!  In put data  files :
!       1.   
!       2.
!  Output data files  :
!       1. Im_result.txt计算的中间数据
!       2.  result.txt计算结果
!-----------------------------------------------------


use jacobi
implicit real*8(a-z)

integer,parameter::N=3

real*8 ::A(N,N),b(N),x(N),x0(N)

  open(unit=101,file='result.txt')
  open(unit=102,file='Im_result.txt')

  x0=(/0d0,0d0,0d0/)

  b=(/20d0,33d0,12d0/)

  A=reshape((/8,4,2,-3,11,1,2,-1,4 /),(/3,3/))

  call solve(A,b,x,x0,N)

  write(101,501)x
  501 format(/,T20,'jacobi迭代法',/,T10,'x(1)',T30,'x(2)',T50,'x(3)',/,3F20.15)

end program main
```

计算结果如下：

$x(1)=2.999999980059588$

$x(2)=2.000000028721297$

$x(3)=1.000000032806938$

迭代过程如下：

| 序号 |   $x(1)$   |   $x(2)$   |   $x(3)$   |
| :--: | :--------: | :--------: | :--------: |
|  1   | 2.50000000 | 3.00000000 | 3.00000000 |
|  2   | 2.87500000 | 2.36363636 | 1.00000000 |
|  3   | 3.13636364 | 2.04545455 | 0.97159091 |
|  4   | 3.02414773 | 1.94783058 | 0.92045455 |
|  5   | 3.00032283 | 1.98398760 | 1.00096849 |
|  6   | 2.99375323 | 1.99997065 | 1.00384168 |
|  7   | 2.99902857 | 2.00262080 | 1.00313072 |
|  8   | 3.00020012 | 2.00063786 | 0.99983051 |
|  9   | 3.00028157 | 1.99991182 | 0.99974048 |
|  10  | 3.00003181 | 1.99987402 | 0.99988126 |
|  11  | 2.99998244 | 1.99997764 | 1.00001559 |
|  12  | 2.99998772 | 2.00000780 | 1.00001437 |
|  13  | 2.99999933 | 2.00000577 | 1.00000419 |
|  14  | 3.00000112 | 2.00000062 | 0.99999889 |
|  15  | 3.00000051 | 1.99999949 | 0.99999929 |
|  16  | 2.99999999 | 1.99999975 | 0.99999987 |
|  17  | 2.99999994 | 1.99999999 | 1.00000007 |

<font face= "times new roman" color=blue > 1.2 Gauss-Seidel迭代法 </font>

由于雅克比迭代法对已经算出来的信息未加充分利用，在计算$x_2$时$x_1$已经计算出来了，计算$x_i$时$x_1,x_2,\cdots ,x_{i-1}$都已经算出，后面的计算值$x_i^{(k+1)}$比前面的计算值$x_i^{(k)}​$要精确些，对<font face= "times new roman" > Jacobi </font>迭代法进行改进：
$$
x_i^{(k+1)}=x_i^{(k)}+\frac{1}{a_{ii}}(b_i-\sum_{j=1}^{i-1}a_{ij}x_j^{(k+1)}-\sum_{j=i}^na_{ij}x_j^{(k)}),(i=1,2,\cdots,n)
\tag{10}
$$

进一步得到其迭代矩阵为：
$$
B_G=(D-L)^{-1}U,f_G=(D-L)^{-1}b\tag{11}
$$

例：用<font face="times new roman">G-S迭代法求解下述方程组</font>
$$
\begin{bmatrix}
8&-3&2\\4&11&-1\\2&1&4
\end{bmatrix}\begin{bmatrix}x_1\\x_2\\x_3\end{bmatrix}=\begin{bmatrix}20\\33\\12\end{bmatrix}\tag{12}
$$
该方程有精确解$X=(3,2,1)$

使用fortran程序求解，代码如下：

代码中使用的迭代方程为公式（9）的变形，如下：
$$
x_i^{(k+1)}=\frac{1}{a_{ii}}(b_i-\sum_{j=1}^{i-1}a_{ij}x_j^{(k+1)}-\sum_{j=i-1}^na_{ij}x_j^{(k)}),(i=1,2,\cdots,n)
\tag{13}
$$

```fortran
module gs
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010-4-5
!-----------------------------------------------------
!  Description : GS迭代法模块
!    
!-----------------------------------------------------
!  Parameters  :
!      1.    IMAX最大允许迭代次数  
!      2.    tol误差容限
!  Contains    :
!      1.    solve GS迭代法方法函数
!      2.
!-----------------------------------------------------
!  Post Script :
!      1.
!      2. 
!-----------------------------------------------------

implicit real*8(a-z)
integer::IMAX=200
real*8::tol=1d-7

contains

subroutine solve(A,b,x,x0,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  GS迭代法函数
!               用于计算方程 AX=b
!-----------------------------------------------------
!  Input  parameters  :
!       1.  A,b 意义即  AX=b
!       2.  x0迭代初值
!       3.  N 方程的维数
!  Output parameters  :
!       1. x 方程的解
!       2.
!  Common parameters  :
!
!----------------------------------------------------

implicit real*8(a-z)
integer::N
integer::i,j,k

real*8::A(N,N),b(N),x(N),x0(N)

real*8::x1(N),x2(N)



!写入标题
  write(102,501)
  501 format(//,18x,'G-S迭代法',//)

!迭代之前两值都设为初值
x1=x0
x2=x1

do k=1,IMAX
   
   do i=1,N
        s=0
        
        do j=1,N
       !-----------------------
       !这段为GS迭代法的核心部分
       !如果j<i 则表示这些量已经更新过了，则下一个元素就用最新的量计算
       !如果j>i 则还没有计算到这些量，所以就用上一次迭代的结果
        if (j<i) then
         s=s+A(i,j)*x2(j)
        else if (j>i) then
        s=s+A(i,j)*x1(j)
        end if
        
        end do 
       !------------------------
       x2(i)=(b(i)-s)/A(i,i)  
   
   end do
   

 !这段程序用于判断精度，满足精度时退出循环,使用向量的二范数继续精度判断.   
   dx2=0
   do i=1,N
    dx2=dx2+(x1(i)-x2(i))**2
   end do
   dx2=dsqrt(dx2)
   
   if (dx2<tol)  exit
!----------------------------------------      
   x1=x2
  
  !记录迭代中间值
     write(102,502)k,x1
     502 format(I3,3F12.8)
  !----
end do 

x=x2

end subroutine solve


end module gs



program  main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-5
!-----------------------------------------------------
!  Purpose   :  采用G-S迭代法计算线性方程
!    
!-----------------------------------------------------
!  In put data  files :
!       1.   
!       2.
!  Output data files  :
!       1. Im_result.txt计算的中间数据
!       2.  result.txt计算结果
!-----------------------------------------------------


use gs
implicit real*8(a-z)

integer,parameter::N=3

real*8 ::A(N,N),b(N),x(N),x0(N)

  open(unit=101,file='result.txt')
  open(unit=102,file='Im_result.txt')

  x0=(/0d0,0d0,0d0/)

  b=(/20d0,33d0,12d0/)

  A=reshape((/8,4,2,-3,11,1,2,-1,4/),(/3,3/))
  
  call solve(A,b,x,x0,N)
  
  write(101,501)x                                                                                                            
  501 format(/,T20,'G-S迭代法',/,T10,'x(1)',T30,'x(2)',T50,'x(3)',/,3F20.15)

end program main
  
```

计算结果如下：

$x(1)=3.000000006322257$

$x(2)=1.999999998008782$

$x(3)=0.999999997336676$

迭代过程如下：

| 序号 |   $x(1)$   |   $x(2)$   |   $x(3)$   |
| :--: | :--------: | :--------: | :--------: |
|  1   | 2.50000000 | 2.09090909 | 1.22727273 |
|  2   | 2.97727273 | 2.02892562 | 1.00413223 |
|  3   | 3.00981405 | 1.99680691 | 0.99589125 |
|  4   | 2.99982978 | 1.99968838 | 1.00016302 |
|  5   | 2.99984239 | 2.00007213 | 1.00006077 |
|  6   | 3.00001186 | 2.00000121 | 0.99999377 |
|  7   | 3.00000201 | 1.99999870 | 0.99999932 |
|  8   | 2.99999968 | 2.00000005 | 1.00000014 |
|  9   | 2.99999998 | 2.00000002 | 1.00000000 |

可以发现，在达到相同精度的情况下，使用<font face="times new roman">G-S</font>迭代法比<font face="times new roman">Jacobi</font>迭代法的迭代次数更少.但在值得指出的是在迭代的收敛性方面<font face="times new roman">G-S</font>迭代法并不比比<font face="times new roman">Jacobi</font>迭代法更具优势，有的时候<font face="times new roman">G-S</font>迭代法可以收敛而<font face="times new roman">Jacobi</font>迭代法不收敛，有的时候<font face="times new roman">Jacobi</font>迭代法可以收敛而<font face="times new roman">G-S</font>迭代法不收敛.

一般地的迭代矩阵，我们可以判断其谱半径是否小于1来判读该迭代方法是否收敛.

> 对于任意$x^{(0)}$和$f$均收敛的充要条件是$\rho (B)<1$.
> 

<font face= "times new roman" color=blue > 1.3 逐次超松弛迭代法（SOR-Successive Over Relaxation） </font>

使用迭代法的困难所在是计算量难以估计. 有时迭代过程虽然收敛，但由于收敛速度缓慢，使计算量变的很大而失去使用价值. 因此，迭代过程的加速具有重要意义. 前面介绍了求解线性方程组的<font face="times new roman">G-S</font>迭代法，通过对其进行<font face="times new roman">SOR</font>加速可得到其他迭代法.通常加速后的迭代法具有更快的收敛速度.

逐次超松弛迭代的迭代公式如下：
$$
x_i^{(k+1)}=(1-\omega)x_i^{(k)}+\frac{\omega}{a_{ii}}(b_i-\sum_{j=1}^{i-1}a_{ij}x_j^{{(k+1)}}-\sum_{j=i+1}^na_{ij}x_j^{(k)}),(j=1,2,\cdots,n)\tag{14}
$$
其中$\omega$为松弛因子，其取值范围$0\leq \omega \leq 1$，当$\omega=1$时，它就是<font face="times new roman">G-S</font>迭代法. 实际上可以证明<font face="times new roman">SOR</font>收敛的必要条件为$0<\omega < 2$ .当松弛因子大于1时称为逐
次超松弛迭代，松弛因子小于1 时称为逐次低松弛迭代。不过，有时都统称为逐次超松弛迭代，或简称<font face="times new roman">SOR</font>迭代方法. 在实际计算中，最优松弛因子很难事先确定，需要通过经验，甚至需要试算才行. 目前已针对某些特殊矩阵建立了最佳松弛因子公式.

进一步得到其迭代矩阵：

$B_\omega=(D-\omega L)^{-1}[(1-\omega)D+\omega U]，f_\omega=\omega(D-\omega L)^{-1}b$

例：用<font face="times new roman">SOR迭代法求解下述方程组</font>
$$
\begin{bmatrix}
-4&1&1&1\\1&-4&1&1\\1&1&-4&1\\1&1&1&-4
\end{bmatrix}\begin{bmatrix}x_1\\x_2\\x_3\\x_4\end{bmatrix}=\begin{bmatrix}1\\1\\1\\1\end{bmatrix}\tag{12}
$$
该方程有精确解$X=(-1,-1,-1,-1)$

使用fortran程序求解，代码如下：

```fortran
module sor
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010-4-5
!-----------------------------------------------------
!  Description : SOR迭代法模块
!    
!-----------------------------------------------------
!  Parameters  :
!      1.    IMAX最大允许迭代次数  
!      2.    tol误差容限
!      3.    omiga 松弛因子
!  Contains    :
!      1.    solve SOR迭代法方法函数
!      2.
!-----------------------------------------------------
!  Post Script :
!      1.
!      2. 
!-----------------------------------------------------

implicit real*8(a-z)
integer::IMAX=200
real*8::tol=1d-5
real*8::omiga=1.5!松弛因子在此修改

contains

subroutine solve(A,b,x,x0,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  SOR迭代法函数
!               用于计算方程 AX=b
!-----------------------------------------------------
!  Input  parameters  :
!       1.  A,b 意义即  AX=b
!       2.  x0迭代初值
!       3.  N 方程的维数
!  Output parameters  :
!       1. x 方程的解
!       2.
!  Common parameters  :
!
!----------------------------------------------------

implicit real*8(a-z)
integer::N
integer::i,j,k

real*8::A(N,N),b(N),x(N),x0(N)

real*8::x1(N),x2(N)



!写入标题
  write(102,501)
  501 format(//,18x,'SOR迭代法',//)

!迭代之前两值都设为初值
x1=x0
x2=x1

do k=1,IMAX
   
   do i=1,N
        s=0
        
        do j=1,N
       !-----------------------
       !如果j<i 则表示这些量已经更新过了，则下一个元素就用最新的量计算
       !如果j>i 则还没有计算到这些量，所以就用上一次迭代的结果
        if (j<i) then
         s=s+A(i,j)*x2(j)
        else if (j>i) then
        s=s+A(i,j)*x1(j)
        end if
        
        end do 
       !------------------------
       x2(i)=(b(i)-s)*omiga/A(i,i)+(1-omiga)*x1(i)
   
   end do
   

 !这段程序用于判断精度，满足精度时退出循环   
   dx2=0
   do i=1,N
    dx2=dx2+(x1(i)-x2(i))**2
   end do
   dx2=dsqrt(dx2)
   
   if (dx2<tol)  exit
!----------------------------------------      
   x1=x2
  
  !记录迭代中间值
     write(102,502)k,x1
     502 format(I3,4F12.8)
  !----
end do 

x=x2

end subroutine solve

end module sor



program  main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-5
!-----------------------------------------------------
!  Purpose   :  采用SOR迭代法计算线性方程
!    
!-----------------------------------------------------
!  In put data  files :
!       1.   
!       2.
!  Output data files  :
!       1. Im_result.txt计算的中间数据
!       2.  result.txt计算结果
!-----------------------------------------------------


use sor
implicit real*8(a-z)

integer,parameter::N=4

real*8 ::A(N,N),b(N),x(N),x0(N)

  open(unit=101,file='result.txt')
  open(unit=102,file='Im_result.txt')


  x0=(/0d0,0d0,0d0,0d0/)

  b=(/1,1,1,1/)

  A=reshape((/-4,1,1,1,&
              1,-4,1,1,&
              1,1,-4,1,&
              1,1,1,-4 /),(/4,4/))
  
  call solve(A,b,x,x0,N)
  
  write(101,501)x                                                                                                            
  501 format(/,T10,'SOR迭代法',//,&
                2x,'x(1)=',F15.8,/,&
                2x,'x(2)=',F15.8,/,&
                2x,'x(3)=',F15.8,/,&
                2x,'x(4)=',F15.8,/)


end program main
 
```

计算结果：

| 松弛因子$\omega$ | 迭代次数 | 松弛因子$\omega$ | 迭代次数 |
| :--------------: | :------: | :--------------: | :------: |
|       1.0        |    22    |       1.5        |    17    |
|       1.1        |    17    |       1.6        |    23    |
|       1.2        |    12    |       1.7        |    33    |
|     **1.3**      |  **11**  |       1.8        |    53    |
|       1.4        |    14    |       1.9        |   109    |

从结果可以看出，当松弛因子为1.3时，迭代次数最少.

其他的迭代方法还有Richardson迭代法，Jacobi超松弛迭代法，最速下降法，共轭梯度法等. 其中Richardson迭代法，Jacobi超松弛迭代法构造方法比较朴素.

SOR迭代法是解大型稀疏线性方程组的有效方法，但需要选取合适的松弛因子，后来发展的共轭梯度法也是求解大型稀疏方程组的理想方法，它比最速下降法的收敛速度快得多，该方法的关键在于寻找共轭方向. 该部分涉及到变分理论和优化方法的思想，以后学会了再更新此部分内容.

参考书目:
[1]张宏伟 计算机科学计算
[2]宋叶志 Fotran95/2003 科学计算与工程
[3]李庆扬 数值分析
		


