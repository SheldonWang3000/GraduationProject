#include <cstdio>
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
using namespace std;
struct pp
{
	int x;
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
	return a.x > b.x;
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
	vector<int> list;
	for (int i = 0; i < 10; ++ i)
	{
		list.push_back(i);
	}
	system("pause");
	return 0;
}