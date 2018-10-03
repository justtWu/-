#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include<cstdio>
//#include <NTL/mat_ZZ_p.h>
#include <NTL/mat_ZZ.h>
//#include <NTL/ZZ_p.h>
#include <NTL/ZZ.h>
//#include <NTL/vec_ZZ_p.h>
#include <NTL/vec_ZZ.h>
#include <NTL/vec_RR.h>
#include <NTL/mat_RR.h>
#include <cmath>
#include <stack>
#include<time.h>
#include <stdlib.h>
#include<math.h>
#include <unistd.h>
using namespace std;
using namespace NTL;

const ZZ w(1ll << 40ll);//左移40位
ZZ aBound(1000), tBound(aBound), eBound(1000);
int l = 100;

const mat_ZZ hCat(const mat_ZZ& A, const mat_ZZ& B);
const mat_ZZ vCat(const mat_ZZ& A, const mat_ZZ& B);

const vec_ZZ decrypt(const mat_ZZ& S, const vec_ZZ& c);

// returns c*
const vec_ZZ getBitVector(const vec_ZZ& c);

// returns S*
const mat_ZZ getBitMatrix(const mat_ZZ& S);

// returns S
const mat_ZZ getSecretKey(const mat_ZZ& T);

// returns M
const mat_ZZ keySwitchMatrix(const mat_ZZ& S, const mat_ZZ& T);

// finds c* then returns Mc*
const vec_ZZ keySwitch(const mat_ZZ& M, const vec_ZZ& c);

// as described, treating I as the secret key and wx as ciphertext
const vec_ZZ encrypt(const mat_ZZ& T, const vec_ZZ& x);

const mat_ZZ getRandomMatrix(long row, long col, const ZZ& bound);

// server side addition with same secret key
const vec_ZZ addn(const vec_ZZ& c1, const vec_ZZ& c2);

// server side linear transformation,
// returns S(Gx) given c=Sx and M (key switch matrix from GS to S)
const vec_ZZ linearTransform(const mat_ZZ& M, const vec_ZZ& c);

// returns M, the key switch matrix from GS to S,
// to be sent to server
const mat_ZZ linearTransformClient(const mat_ZZ& T, const mat_ZZ& G);

// computes an inner product, given two ciphertexts and the keyswitch matrix
const vec_ZZ innerProd(const vec_ZZ& c1, const vec_ZZ& c2, const mat_ZZ& M);

// returns M, the key switch matrix from vec(S^t S) to S,
// to be sent to the server
const mat_ZZ innerProdClient(const mat_ZZ& T);

// returns a column vector
const mat_ZZ vectorize(const mat_ZZ& M);

const mat_ZZ copyRows(const mat_ZZ& row, long numrows);




// finds c* then returns Mc*
const vec_ZZ keySwitch(const mat_ZZ& M, const vec_ZZ& c) {
	vec_ZZ cstar = getBitVector(c);
	return M * cstar;
}


const mat_ZZ getRandomMatrix(long row, long col, const ZZ& bound){
	mat_ZZ A;
	A.SetDims(row, col);
	for (int i=0; i<row; ++i){
		for (int j=0; j<col; ++j){
			A[i][j] = RandomBnd(bound);
		}
	}
	return A;
}




// returns S*
const mat_ZZ getBitMatrix(const mat_ZZ& S) {
	mat_ZZ result;
	int rows = S.NumRows(), cols = S.NumCols();
	result.SetDims(rows, l * cols);

	vec_ZZ powers;
	powers.SetLength(l);
	powers[0] = 1;
	for(int i = 0; i < l - 1; ++i) {
		powers[i+1] = powers[i]*2;//1-2**99
	}

	for(int i = 0; i < rows; ++i) {
		for(int j = 0; j < cols; ++j) {
			for(int k = 0; k < l; ++k) {
				result[i][j*l + k] = S[i][j] * powers[k];
			}
		}
	}

	return result;
}


// returns c*
const vec_ZZ getBitVector(const vec_ZZ& c) {
	vec_ZZ result;
	int length = c.length();//length=3472
	result.SetLength(length * l);
	for(int i = 0; i < length; ++i) {
		ZZ sign = (c[i] < ZZ(0)) ? ZZ(-1) : ZZ(1);
		ZZ value = c[i] * sign;//值都变为正数
		for(int j = 0; j < l; ++j) {
			result[i * l + j] = sign*bit(value, j);
		}
	}
	return result;//长度为347200的比特串
}



// returns S,扩充矩阵T,行数不变，列数等于I+T
const mat_ZZ getSecretKey(const mat_ZZ& T) {
	mat_ZZ I;
	ident(I, T.NumRows());//生成单位矩阵I，I的行数为T的行数
	return hCat(I, T);//以列合并I与T
}

//合并矩阵A，B
const mat_ZZ hCat(const mat_ZZ& A, const mat_ZZ& B) {
	assert(A.NumRows() == B.NumRows());

	int rows = A.NumRows(), colsA = A.NumCols(), colsB = B.NumCols();
	mat_ZZ result;
	result.SetDims(rows, colsA + colsB);

	// Copy A
	for(int i = 0; i < rows; ++i) {
		for(int j = 0; j < colsA; ++j) {
			result[i][j] = A[i][j];
		}
	}

	// Copy B
	for(int i = 0; i < rows; ++i) {
		for(int j = 0; j < colsB; ++j) {
			result[i][colsA + j] = B[i][j];
		}
	}

	return result;
}

//列数不变，以行合并A和B
const mat_ZZ vCat(const mat_ZZ& A, const mat_ZZ& B) {
	assert(A.NumCols() == B.NumCols());

	int cols = A.NumCols(), rowsA = A.NumRows(), rowsB = B.NumRows();
	mat_ZZ result;
	result.SetDims(rowsA + rowsB, cols);

	// Copy A
	for(int i = 0; i < rowsA; ++i) {
		for(int j = 0; j < cols; ++j) {
			result[i][j] = A[i][j];
		}
	}

	// Copy B
	for(int i = 0; i < rowsB; ++i) {
		for(int j = 0; j < cols; ++j) {
			result[i + rowsA][j] = B[i][j];
		}
	}

	return result;
}
//返回最近的整数
inline const ZZ nearestInteger(const ZZ& x, const ZZ& w) {
	return (x + (w+1)/2) / w;
}

const vec_ZZ decrypt(const mat_ZZ& S, const vec_ZZ& c) {
	vec_ZZ Sc = S*c;
	vec_ZZ output;
	output.SetLength(Sc.length());
	for (int i=0; i<Sc.length(); i++) {
		output[i] = nearestInteger(Sc[i], w);
	}
	return output;
}

const mat_ZZ keySwitchMatrix(const mat_ZZ& S, const mat_ZZ& T) {
	mat_ZZ Sstar = getBitMatrix(S);//列数扩充100倍，并且数字以2的指数的倍数扩大
	mat_ZZ A = getRandomMatrix(T.NumCols(),Sstar.NumCols(),aBound);//随机操作
	mat_ZZ E = getRandomMatrix(Sstar.NumRows(),Sstar.NumCols(),eBound);//随机操作
	return vCat(Sstar + E - T*A, A);//以行合并参数1和参数2
}

const vec_ZZ encrypt(const mat_ZZ& T, const vec_ZZ& x) {
	mat_ZZ I;
	ident(I, x.length());
	return keySwitch(keySwitchMatrix(I, T), w * x);
}




const vec_ZZ addVectors(const vec_ZZ& c1, const vec_ZZ& c2){
	return c1 + c2;
}

const vec_ZZ linearTransform(const mat_ZZ& M, const vec_ZZ& c){
	return M * getBitVector(c);
}

const mat_ZZ linearTransformClient(const mat_ZZ& G, const mat_ZZ& S, const mat_ZZ& T){
	return keySwitchMatrix(G * S, T);
}


const vec_ZZ innerProd(const vec_ZZ& c1, const vec_ZZ& c2, const mat_ZZ& M){
	mat_ZZ cc1;
	mat_ZZ cc2;
	mat_ZZ cc;

	cc1.SetDims(c1.length(), 1);
	for (int i=0; i<c1.length(); ++i){
		cc1[i][0] = c1[i];
	}
	cc2.SetDims(1, c2.length());
	for (int i=0; i<c2.length(); ++i){
		cc2[0][i] = c2[i];
	}
	cc = vectorize(cc1 * cc2);

	vec_ZZ output;
	output.SetLength(cc.NumRows());
	for (int i=0; i<cc.NumRows(); i++) {
		output[i] = nearestInteger(cc[i][0], w);
	}
	return M * getBitVector(output);
}

const mat_ZZ innerProdClient(const mat_ZZ& T,const mat_ZZ& S){
	//mat_ZZ S = getSecretKey(T);
	mat_ZZ tvsts = transpose(vectorize(transpose(S) * S));
	mat_ZZ mvsts = copyRows(tvsts, T.NumRows());
	return keySwitchMatrix(mvsts, T);
}




const mat_ZZ copyRows(const mat_ZZ& row, long numrows){
	mat_ZZ ans;
	ans.SetDims(numrows, row.NumCols());
	for (int i=0; i<ans.NumRows(); ++i){
		for (int j=0; j<ans.NumCols(); ++j){
			ans[i][j] = row[0][j];
		}
	}
	return ans;
}

const mat_ZZ vectorize(const mat_ZZ& M){
	mat_ZZ ans;
	ans.SetDims(M.NumRows() * M.NumCols(), 1);
	for (int i=0; i<M.NumRows(); ++i){
		for (int j=0; j<M.NumCols(); ++j){
			ans[i*M.NumCols() + j][0] = M[i][j];
		}
	}
	return ans;
}

//************关于fisherLMC自定义函数*********

struct data
{
	mat_RR dataSet;
	mat_RR dataTag;
	int rows;
};
//the fisherLMC for cleartext with division

class flmcOnClearText
{
	private:
		data trainData;	
		mat_RR w0;
		mat_RR w_star;
	public:
		flmcOnClearText()
		{
			w_star.SetDims(4,1);
		}

		void trainLMC(data theData)
		{
			trainData=theData;
			mat_RR dataPositive;
			mat_RR dataMinus;
			dataPositive.SetDims(1,4);
			dataMinus.SetDims(1,4);
			mat_RR cPositive;
			mat_RR cMinus;
			mat_RR c;
			cPositive.SetDims(4,4);
			cMinus.SetDims(4,4);
			c.SetDims(4,4);
			RR positiveCount;
			RR minusCount;
			positiveCount = RR(0);
			minusCount = RR(0);

			for(int i=0;i<1097;i++)
			{
				if (trainData.dataTag[i][0]==0)
				{
					dataPositive[0] = dataPositive[0] + trainData.dataSet[i];
					positiveCount++;
				}
				else
				{
					dataMinus[0] = dataMinus[0] + trainData.dataSet[i];
					minusCount++;

				}
			}

			for(int i=0;i<4;i++)
			{
				dataPositive[0][i] = dataPositive[0][i]/positiveCount;
			}

			for(int i=0;i<4;i++)
			{
				dataMinus[0][i] = dataMinus[0][i]/minusCount;
			}
			mat_RR tmp;
			tmp.SetDims(1,4);
			for(int i=0;i<trainData.rows;i++)
			{
				if (trainData.dataTag[i][0]==2)
				{
					tmp[0] = trainData.dataSet[i];
					cPositive = cPositive + (transpose(tmp) - transpose(dataPositive)) * (tmp - dataPositive);
				}
				else
				{
					tmp[0] = trainData.dataSet[i];
					cMinus = cMinus + (transpose(tmp) - transpose(dataMinus)) * (tmp - dataMinus);
				}
				
			}
			c= cPositive + cMinus;
			mat_RR cInverse;
			cInverse = inv(c);
			w_star = cInverse*(transpose(dataPositive) - transpose(dataMinus));
			w0 = -0.5*(transpose(w_star)*transpose(dataPositive) - transpose(w_star)*transpose(dataMinus));
		}

		void testLMC(data testData)
		{
			mat_RR f;
			int errCount=0;
			mat_RR x;
			x.SetDims(1,4);
			int flag=0;
			for(int i=0;i<testData.rows;i++)
			{
				x[0]=testData.dataSet[i];
				f =  transpose(w_star)*transpose(x)- w0;
				//cout<<f<<endl;
				if (f[0][0]<0)
				{
					flag=1;
				}
				else
				{
					flag=0;
				}
				if (flag!=testData.dataTag[i][0])
				{
					errCount++;
				}
			}
			float errRate;
			errRate=float(errCount)/float(testData.rows);
			cout<<"cleartext errRate: "<<errRate<<endl;
		}

};



int main()
{
	data trainData;
	data testData;

	trainData.dataSet.SetDims(1097,4);
	trainData.dataTag.SetDims(1097,1);
	trainData.rows=1097;

	testData.dataSet.SetDims(275,4);
	testData.dataTag.SetDims(275,1);
	testData.rows=275;



	mat_RR traindata;//训练数据集 
	traindata.SetDims(1097,4);
	mat_RR testdata;
	testdata.SetDims(275,4);


	vec_RR trainTag;//读入训练数据集的标签 
	trainTag.SetLength(1097);
	vec_RR testTag;//读入训练数据集的标签 
	testTag.SetLength(275);



	freopen("D_train.txt", "r", stdin);
	for(int i=0; i<1097;i++)
	{
		for(int j=0; j<4; j++)
		{
			cin>>trainData.dataSet[i][j];		
		}
		
	}
	//cout<<trainData.dataSet<<endl;


	freopen("D_train_tag.txt", "r", stdin);
	for(int i=0; i<1097;i++)
	{
		cin>>trainData.dataTag[i][0];	
	}
	//cout<<trainData.dataTag<<endl;
	freopen("D_test.txt", "r", stdin);
	for(int i=0; i<275;i++)
	{
		for(int j=0; j<4; j++)
		{
			cin>>testData.dataSet[i][j];
		}
		
	}
	//cout<<testData.dataSet<<endl;

	fstream fp_S("D_test_tag.txt",fstream::in) ;
	assert(fp_S.is_open());
	for(int i=0;i<275;i++)
	{	
		fp_S >> testData.dataTag[i][0];	
	}
	fp_S.close();

	flmcOnClearText fLMC;
	fLMC.trainLMC(trainData);
	fLMC.testLMC(testData);
}