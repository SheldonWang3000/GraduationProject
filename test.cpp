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
	const int MAXN = 20;
	int a[MAXN];
	int i;
	for (i = 0; i < MAXN; ++i)
		a[i] = rand() % (MAXN * 2);

	//��̬����vector ����vector����
	vector<pp> *pvet = new vector<pp>(20);
	//pvet->assign(a, a + MAXN);
	for (int i = 0; i < MAXN; ++ i)
	{
		(*pvet)[i].x = a[i];
	}
	
	PrintfVectorInt(*pvet);
	//����
	make_heap(pvet->begin(), pvet->end(), compare);
	PrintfVectorInt(*pvet);

	pop_heap(pvet->begin(), pvet->end(), compare);
	pop_heap(pvet->begin(), -- pvet->end(), compare);
	PrintfVectorInt(*pvet);
	int temp = (*pvet)[pvet->size() - 1].x;
	pvet->pop_back();
	temp += (*pvet)[pvet->size() - 1].x;
	pvet->pop_back();

	pp kk;
	kk.x = temp;
	pvet->push_back(kk);
	PrintfVectorInt(*pvet);
	push_heap(pvet->begin(), pvet->end(), compare);
	PrintfVectorInt(*pvet);
	//PrintfVectorInt(*pvet);

	////���������� ���������м��룬�ٵ���push_heap()
	//pvet->push_back(25);
	//push_heap(pvet->begin(), pvet->end());
	//PrintfVectorInt(*pvet);

	////ɾ������  Ҫ�ȵ���pop_heap()������������ɾ��
	//pop_heap(pvet->begin(), pvet->end());
	//pvet->pop_back();
	//pop_heap(pvet->begin(), pvet->end());
	//pvet->pop_back();
	//PrintfVectorInt(*pvet);

	////������
	//sort_heap(pvet->begin(), pvet->end());
	//PrintfVectorInt(*pvet);

	delete pvet;
	system("pause");
	return 0;
}