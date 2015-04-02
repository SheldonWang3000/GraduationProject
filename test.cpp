#include <cstdio>
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <time.h>
#include "point_data.h"
using namespace std;
struct pp
{
	int x;
	friend bool operator==(const pp &l, const pp &r)
	{
		return l.x == r.x;
	}
};

struct link
{
	int a;
	link *next;
};
void PrintfVectorInt(vector<pp> &vet)
{
	for (vector<pp>::iterator pos = vet.begin(); pos != vet.end(); pos++)
		printf("%d ", pos->x);
	putchar('\n');
}

bool compare(pp a, pp b)
{
	//a < b ������, a > b ����С��
	return a.x < b.x;
}
int main()
{
	//const int MAXN = 20;
	//int a[MAXN];
	//int i;
	//for (i = 0; i < MAXN; ++i)
	//	a[i] = rand() % (MAXN * 2);

	////��̬����vector ����vector����
	//vector<pp> *pvet = new vector<pp>(20);
	////pvet->assign(a, a + MAXN);
	//for (int i = 0; i < MAXN; ++ i)
	//{
	//	(*pvet)[i].x = a[i];
	//}
	//
	//PrintfVectorInt(*pvet);
	////����
	//make_heap(pvet->begin(), pvet->end(), compare);
	//PrintfVectorInt(*pvet);

	//pop_heap(pvet->begin(), pvet->end(), compare);
	//pop_heap(pvet->begin(), -- pvet->end(), compare);
	//PrintfVectorInt(*pvet);
	//int temp = (*pvet)[pvet->size() - 1].x;
	//pvet->pop_back();
	//temp += (*pvet)[pvet->size() - 1].x;
	//pvet->pop_back();

	//pp kk;
	//kk.x = temp;
	//pvet->push_back(kk);
	//PrintfVectorInt(*pvet);
	//push_heap(pvet->begin(), pvet->end(), compare);
	//PrintfVectorInt(*pvet);
	////PrintfVectorInt(*pvet);

	//////���������� ���������м��룬�ٵ���push_heap()
	////pvet->push_back(25);
	////push_heap(pvet->begin(), pvet->end());
	////PrintfVectorInt(*pvet);

	//////ɾ������  Ҫ�ȵ���pop_heap()������������ɾ��
	////pop_heap(pvet->begin(), pvet->end());
	////pvet->pop_back();
	////pop_heap(pvet->begin(), pvet->end());
	////pvet->pop_back();
	////PrintfVectorInt(*pvet);

	//////������
	////sort_heap(pvet->begin(), pvet->end());
	////PrintfVectorInt(*pvet);

	//delete pvet;
	
	//int a[10];
	//for (int i = 0; i != 10; ++ i)
	//{
	//	a[i] = i;
	//}
	//link *head = nullptr;
	//link *last;
	//link *temp;
	//for (int i = 0; i < 10; ++ i)
	//{
	//	
	//	if (head == nullptr)
	//	{
	//		head = new link;
	//		head->a = i;
	//		head->next = nullptr;
	//		last = head;
	//	}
	//	else
	//	{
	//		temp = new link;
	//		temp->a = i;
	//		temp->next = nullptr;
	//		last->next = temp;
	//		last = last->next;
	//	}
	//}
	//while (head->next)
	//{
	//	cout << head->a << endl;
	//	head = head->next;
	//}
	//cout << head->a << endl;
	vector<pp> list;
	for (int i = 0; i < 10000; ++i)
	{
		pp temp;
		temp.x = i;
		list.push_back(temp);
	}
	make_heap(list.begin(), list.end(), compare);
	clock_t start = clock();
	for (int i = 10000; i < 12000; ++i)
	{
		pp temp;
		temp.x = i;
		list.push_back(temp);
	}
	make_heap(list.begin(), list.end(), compare);
	clock_t end = clock();
	cout << list[0].x << endl;
	cout << (double)(end - start) / CLOCKS_PER_SEC << endl;

	list.clear();
	for (int i = 0; i < 10000; ++i)
	{
		pp temp;
		temp.x = i;
		list.push_back(temp);
	}
	make_heap(list.begin(), list.end(), compare);
	
	start = clock();
	for (int i = 10000; i < 12000; ++i)
	{
		pp temp;
		temp.x = i;
		list.push_back(temp);
		push_heap(list.begin(), list.end(), compare);
	}
	end = clock();
	cout << list[0].x << endl;
	cout << (double)(end - start) / CLOCKS_PER_SEC << endl;
	system("pause");
	return 0;
}