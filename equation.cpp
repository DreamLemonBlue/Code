#include<iostream.h>
#include<math.h>
#include<stdlib.h>
#include<fstream.h>																						 //函数头


/*-------------------------------------------------------------------------------------*/
//在这里修改你的精度

double eps=1e-6;

/*-------------------------------------------------------------------------------------*/
//在这里修改你的函数f(x)

double f(double x)
{
	return (4 * powl(x, 3) - 8 * powl(x, 2) + 2 * x + 1);
}

/*-------------------------------------------------------------------------------------*/
//在这里修改你的导数df(x)

double df(double x)
{
	return (12 * powl(x, 2) - 16 * x + 2);
}

/*-------------------------------------------------------------------------------------*/
//在这里修改你的迭代式g(x)

double g(double x)
{
	double p=1.0/3;
	return (powl((8*powl(x,2)-2*x-1),p));
}

/*-------------------------------------------------------------------------------------*/
//以下代码请勿修改
//函数声明区

void welcome();																								//欢迎页面

class solve;																										//声明类

void method_select(solve &a);																	//方法分类

double half(double a,double b,double eps);												//二分法算法函数声明
double newn(double x0, double eps);															//牛顿迭代法算法函数声明
double itr(double (*fun)(double), double x0, double eps);						//简单迭代法算法函数声明
double chord(double x0, double eps);														//弦截法算法函数声明
double stef(double x0, double x);																//史蒂芬孙算法组件函数声明
double steffensen(double x0, double eps);												//史蒂芬孙算法函数声明

double gc(double x0,double x)																		//弦截法公式
{
	return (x-f(x)/(f(x)-f(x0))*(x-x0));
}
double g1(double x)																						//简单迭代公式1
{
	return (-1.0) / (4 * pow(x, 2) - 8 * x + 2);
}
double g2(double x)																						//简单迭代公式2
{
	double t=1.0/2;
	return pow(((4*pow(x,3)+2*x+1)*(1.0/8)),t);
}
double g3(double x)																						//简单迭代公式3,即上面的迭代式g(x)
{
	double t=1.0/3;
	return pow(((8*pow(x,2)-2*x-1)*(1.0/4)),t);
}
double g4(double x0)																					//牛顿迭代法逼近
{
	return (x0 - f(x0) / df(x0));
}
double g5(double x)																						//史蒂芬孙算法迭代公式
{
	return ((-4 * pow(x, 3) + 8 * pow(x, 2) - 1) * (1.0 / 2));
}

class solve																										//解法类
{
	double a,b;																									//两个区间端点
	double h;																									//步长值
	double c[5][2];																							//区间端点储存
	int i;																												//计数变量
	double result[3];																						//结果储存
public:
	solve(double p=0,double q=1,double r=0.1)												//构造函数
	{
		a=p;
		b=q;
		h=r;
	}
	void filter()																									//筛选区间并储存
	{
		i=0;
		double a0,b0;
		a0=a;
		b0=a+h;
		do
		{
			if(f(a0)*f(b0)<0)
			{
				c[i][0]=a0;
				c[i][1]=b0;
				i++;
			}
			a0=b0;
			b0+=h;	
		}while(b0<=b);
	}
	void halfs()																									//二分法
	{
		filter();
		double x;
		cout<<"IT\ta1\tb1\tx\tf(x)\n";
		for(int j=0;j<i;j++)
		{
			a=c[j][0];
			b=c[j][1];
			x=half(a,b,eps);
			result[j] = x;
			cout<<"x="<<x<<endl;
		}
	}
	void itrs()																									//简单迭代法
	{
		filter();
		double x;
		cout << "IT\ta1\tb1\tx\tf(x)\n";
		x = itr(g1, c[0][0], eps);
		result[0] = x;
		cout << "x1=" << x << endl;
		x = itr(g2, c[1][0], eps);
		result[1] = x;
		cout << "x2=" << x << endl;
		x = itr(g3, c[2][0], eps);
		cout << "x3=" << x << endl;
		result[2]=x;
	}
	void newns()																								//牛顿迭代法
	{
		filter();
		double x;
		for (int j = 0; j < i; j++)
		{
			cout << "IT\tx0\t\tx\n";
			x = newn(c[j][0], eps);
			result[j] = x;
			cout << "\nx=" << x << endl;
		}
	}
	void chords()																								//弦截法
	{
		filter();
		double x;
		cout << "IT\tx0\t\tx\n";
		for (int q = 0; q < i; q++)
		{
			x = chord(c[q][0], eps);
			result[q] = x;
			cout << "\nx=" << x << endl;
		}
	}
	void steffensens()																						//Steffensen法
	{
		filter();
		double x;
		cout << "IT\tx0\t\tx\n";
		for (int r = 0; r < i; r++)
		{
			x = steffensen(c[r][0], eps);
			result[r] = x;
			cout << "x=" << x << endl;
		}
	}
	void output()																								//输出至文件
	{
		ofstream outfile("answer.txt");
		for(int i=0;i<3;i++)
		{
			outfile<<"x"<<i<<"="<<result[i]<<'\n';
		}
		outfile.close();
	}
	
};

void main()																										//主函数
{
	welcome();
	double m,n;																								//区间端点
	cout<<"请输入区间[a,b]端点:"<<endl;
	cout<<"a=";
	cin>>m;
	cout<<"b=";
	cin>>n;
	solve a(m,n);
	method_select(a);
	a.output();
}

void method_select(solve &a)
{
	system("cls");
	cout<<"请输入你想选择方法对应的数字:"<<endl;
	cout<<"1:二分法\n2:简单迭代法\n3:Steffensen法\n4:牛顿法\n5:弦截法\n"<<endl;
	int n;
	cin>>n;
	switch(n)//分支
	{
	case 1:		a.halfs();
		break;
	case 2:		a.itrs();
		break;
	case 3:		a.steffensens();
		break;
	case 4:		a.newns();
		break;
	case 5:		a.chords();
		break;
	}
}

void welcome()
{
	cout<<"欢迎你使用本程序!\n";
	cout<<"本程序由王泽榕，于豪林，冯煜轩共同完成\n";
	cout<<"修改方程请在源文件中修改函数，导函数和迭代函数\n";
	cout<<"默认方程:4x^3-8x^2+2x+1=0\n";
	cout<<"选择方法来解上述方程,得到的解会被保存在answer.txt中\n";
}
//以下是不同算法各自的函数
double half(double a,double b,double eps)
{
	int it=0;
	double x,a1,b1,y0,y;
	a1=a;
	b1=b;
	y0=f(a);
	do
	{
		x=(a1+b1)/2;
		cout<<it<<'\t'<<a1<<'\t'<<b1<<'\t'<<x<<'\t'<<f(x)<<endl;
		y=f(x);
		if(y*y0>0)
		{
			a1=x;
			y0=f(a1);
		}
		else
		{
			b1=x;		
		}
		it++;
	}
	while(fabs(y)>eps);
	return(x);
}//二分法
double newn(double x0, double eps)//牛顿迭代法
{
	int it = 0, e;
	double x;
	do
	{
		x = g4(x0);
		cout << it << '\t' << x0 << '\t' << '\t' << x << endl;
		e = fabs((x - x0) / x) > eps;
		if (e) { x0 = x; it++; }
	} while (e);
	return(x);
}
double itr(double (*fun)(double), double x0, double eps)//简单迭代法
{
	int it = 0, e;
	double x;
	do
	{
		x = fun(x0);
		cout << it << '\t' << x0 << '\t' << '\t' << x << endl;
		e = fabs((x - x0) / x) > eps;
		if (e) { x0 = x; it++; }
	} while (e);
	return(x);
}
double chord(double x0, double eps)//弦截法
{
	int it = 0, e;
	double x, x1;
	if (f(x0) > 0) x1 = x0 + 0.1;
	else x1 = x0 - 0.1;
	do
	{
		x = gc(x0, x1);
		cout << it << '\t' << x0 << '\t' << '\t' << x << endl;
		e = fabs((x1 - x0) / x1) > eps;
		if (e) { x0 = x1; x1 = x; it++; }
	} while (e);
	return(x);
}
double stef(double x0, double x)//Steffensen算法组件
{
	double y, z;
	y = g5(x0);
	z = g5(y);
	return (x0 - (powl((y - x0), 2)) / (z - 2 * y + x0));
}
double steffensen(double x0, double eps)//Stenffensen法
{
	int it = 0, e;
	double x, x1;
	if (f(x0) > 0)
		x1 = x0 + 0.1;
	else
		x1 = x0 - 0.1;
	do
	{
		x = stef(x0, x1);
		cout << it << '\t' << x0 << '\t' << '\t' << x << endl;
		e = fabs((x1 - x0) / x1) > eps;
		if (e)
		{
			x0 = x1;
			x1 = x;
			it++;
		}
	} while (e);
	return(x);
}