#include <stdio.h>

typedef struct
{
	float x;
	float y;
} point;

int main()
{
	point P1, P2;
	scanf("%f%f", &P1.x, &P1.y);
	scanf("%f%f", &P2.x, &P2.y);
	float dx1, dy1;
	dy1 = P1.y - P2.y;
	dx1 = P1.x - P2.x;
	float minx = P2.x;
	float miny = P2.y;
	if(dx1 < 0)
	{
		minx = P1.x;
		dx1 *= -1;
	}
	if(dy1 < 0)
	{
		miny = P1.y;
		dy1 *= -1;
	}
	if(dx1 > dy1)
	{
		for(float i=minx; i<=P1.x || i<=P2.x ; ++i)
		{
			float x = i;
			float y = P1.y + ((x-P1.x)*((float)(P2.y-P1.y))/(P2.x-P1.x));
			int y1 = y;
			if(y > y1 + 0.5)
			{
				y1 = y+1;
			}
			else if(y == y1 + 0.5 && y1%2 !=0)
			{
				y1 = y+1;
			}
			int x1 = x;
			if(x > x1 + 0.5 || (x==x1+0.5 && x1%2!=0))
			{
				x1 = x+1;
			}
			printf("%.2f %.2f\t%d %d\n", x, y, x1, y1);
		}
	}
	else
	{
		for(float i=miny; i<=P1.y || i<=P2.y; ++i)
		{
			float y = i;
			float x = P1.x + ((y-P1.y)*((float)(P2.x-P1.x))/(P2.y-P1.y));
			int x1 = x;
			if(x > x1 + 0.5)
			{
				x1 = x+1;
			}
			else if(x==x1+0.5 && x1%2!=0)
			{
				x1 = x+1;
			}
			int y1 = y;
			if(y > y1 + 0.5 || (y==y1+0.5 && y1%2!=0))
			{
				y1 = y+1;
			}
			printf("%.2f %.2f\t%d %d\n", x, y, x1, y1);
		}
	}
	return 0;
}
