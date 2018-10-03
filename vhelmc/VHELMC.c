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
////////
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

//************关于LMC自定义函数*********

struct data
{
	mat_ZZ dataSet;
	mat_ZZ dataTag;
	int rows;
};
//读数据文件
const data loadDataSet(char* filename,int rows)
{
	//cout<<filename;
	data theData;
	theData.dataSet.SetDims(rows,2);
	theData.dataTag.SetDims(rows,1);

	mat_RR theDataSet;
	mat_RR theDataTag;
	theDataSet.SetDims(rows,2);
	theDataTag.SetDims(rows,1);
	//freopen(filename, "r", stdin);
	fstream fp_S(filename,fstream::in) ;
	assert(fp_S.is_open());
	for(int i=0; i<rows;i++)
	{
		for(int j=0; j<3; j++)
		{
			if(j<2)
			{
				fp_S>>theDataSet[i][j];
				theDataSet[i][j]=theDataSet[i][j]*pow(10,7);
				theData.dataSet[i][j]=TruncToZZ(theDataSet[i][j]);
			}
			else
			{
				fp_S>>theDataTag[i][0];
				theData.dataTag[i][0]=TruncToZZ(theDataTag[i][0]);
			}

			
		}
	}
	fp_S.close();
	//fclose(stdin);
	theData.rows=rows;
	//cout <<"dataSet=" <<theData.dataSet <<endl;
	//cout <<"dataTag=" <<theData.dataTag;
	return theData;
	
}	

//the lmc for cleartext
class lmcOnClearText
{
	private:
		data trainData;
		vec_ZZ dataPositive;
		vec_ZZ dataMinus;	
		int positiveCount;
		int minusCount;
		vec_ZZ x0;
		vec_ZZ y0;
		vec_ZZ theW;	
		ZZ f;
	public:
		lmcOnClearText(data theData)
		{
			dataPositive.SetLength(30);
			dataMinus.SetLength(30);

			x0.SetLength(30);
			theW.SetLength(30);
			
			trainData=theData;

			for(int i=0;i<trainData.rows;i++)
			{
				if (trainData.dataTag[i][0]==2)
				{

					dataPositive = dataPositive + trainData.dataSet[i];
					positiveCount++;
				}
				else
				{
					dataMinus = dataMinus + trainData.dataSet[i];
					minusCount++;
				}
			}
			x0 = minusCount*dataPositive + positiveCount*dataMinus;
			y0 = minusCount*dataPositive - positiveCount*dataMinus;
			dataPositive = dataPositive * minusCount;
			dataMinus = dataMinus * positiveCount;
			theW = ZZ(2)*minusCount*positiveCount*(dataPositive - dataMinus);

		}

		void testLMC(data testData)
		{
			ZZ c;
			c=x0*y0;
			cout<<"c:"<<c<<endl;
			int errCount=0;
			vec_ZZ x;
			x.SetLength(30);
			int flag=0;
			for(int i=0;i<testData.rows;i++)
			{
				x=testData.dataSet[i];
				f =  theW*x- c;
				if (f<0)
				{
					flag=4;
				}
				else
				{
					flag=2;
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

//the lmc for ciphertext
class lmcOnCiphertext
{
	private:
		data trainData;
		vec_ZZ dataPositive;
		vec_ZZ dataMinus;	
		int positiveCount;
		int minusCount;
		vec_ZZ x0;
		vec_ZZ y0;
		vec_ZZ theW;	
		ZZ f;
		mat_ZZ T,S;
	public:
		lmcOnCiphertext(data theData)
		{
			dataPositive.SetLength(3);//after encrypt the length add 1
			dataMinus.SetLength(3);

			x0.SetLength(3);
			y0.SetLength(3);
			theW.SetLength(3);
			
			clear(dataPositive);//初始一定为0 
			clear(dataMinus);

			trainData=theData;
			T = getRandomMatrix(2 , 1 , aBound);//key 
			S = getSecretKey(T);
			vec_ZZ  enc_traindata[90],dec_traindata[90];
			 for(int i=0;i<trainData.rows;i++)
    		{ 
				enc_traindata[i] = encrypt(T, trainData.dataSet[i]);

			}

			for(int i=0;i<trainData.rows;i++)
			{
				if(trainData.dataTag[i][0]==1)
				{
					dataPositive=dataPositive+enc_traindata[i];
					positiveCount++;
				}
				else
				{
					dataMinus=dataMinus+enc_traindata[i];
					minusCount++;
				}
			}
			x0 = minusCount*dataPositive + positiveCount*dataMinus;
			y0 = minusCount*dataPositive - positiveCount*dataMinus;
			dataPositive = dataPositive * minusCount;
			dataMinus = dataMinus * positiveCount;
			theW = ZZ(2)*minusCount*positiveCount*(dataPositive - dataMinus);

		}

		void testLMC(data testData)
		{
			//encrypt test data
			vec_ZZ encryptedData[100];
			for(int i=0;i<testData.rows;i++)
			{
				encryptedData[i]=encrypt(T, testData.dataSet[i]);
			}
			mat_ZZ M0;
			M0 = innerProdClient(T,S);
			vec_ZZ c0;
			vec_ZZ c;
			ZZ c_star;
			cout<<"x*y:"<<x0*y0<<endl;
			c0 = innerProd(y0,x0,M0);//encrypted data innerMutiply
			c_star=c0[0];
			cout<<"c_star:"<<c_star<<endl;
			vec_ZZ x;
			x.SetLength(3);
			int flag=0;
			vec_ZZ enc_tag;
			enc_tag.SetLength(100);

			vec_ZZ temp_enc_tag;
			temp_enc_tag.SetLength(3);

			vec_ZZ temp_dec_tag;
			temp_dec_tag.SetLength(2);

			for(int i=0;i<testData.rows;i++)
			{
				x=encryptedData[i];
				f = theW*x - c_star;
				enc_tag[i] = f; 		
			}
			//cout<<enc_tag<<endl;
			int count=0;
			int tmp=0;
			int errCount=0;
			while(count!=100)
			{
				for(int i=0;i<3;i++)
				{
					if(i<2)
					{

						temp_enc_tag[i]=enc_tag[count];
						count++;
					}
					else
					{
						temp_enc_tag[i]=ZZ(0);
					}
					
				}
				//cout<<temp_enc_tag<<endl;
				temp_dec_tag=decrypt(S,temp_enc_tag);
				//cout<<temp_dec_tag<<endl;
				for(int j=0;j<2;j++)
				{
					if(temp_dec_tag[j]<ZZ(0))
					{
						flag=-1;
					}
					else
					{
						flag=1;
					}
					tmp=count-2+j;
					if(flag!=testData.dataTag[tmp][0])
					{
						errCount++;
					}
					//cout<<flag<<":"<<testData.dataTag[tmp][0]<<endl;
				}	
			}

			float errRate;
			errRate=float(errCount)/float(testData.rows);
			cout<<"ciphertext errRate: "<<errRate<<endl;
		}
		


};


class lmcOnCiphertext2
{
	private:
		data trainData;
		vec_ZZ dataPositive;
		vec_ZZ dataMinus;	
		int positiveCount;
		int minusCount;
		vec_ZZ x0;
		vec_ZZ y0;
		vec_ZZ theW;	
		ZZ f;
		mat_ZZ T,S;
	public:
		lmcOnCiphertext2(data theData)
		{
			dataPositive.SetLength(31);//after encrypt the length add 1
			dataMinus.SetLength(31);

			x0.SetLength(31);
			y0.SetLength(31);
			theW.SetLength(31);
			
			clear(dataPositive);//初始一定为0 
			clear(dataMinus);

			trainData=theData;
			T = getRandomMatrix(30 , 1 , aBound);//key 
			S = getSecretKey(T);
			vec_ZZ  enc_traindata[469],dec_traindata[469];
			 for(int i=0;i<trainData.rows;i++)
    		{ 
				enc_traindata[i] = encrypt(T, trainData.dataSet[i]);
	
			}

			for(int i=0;i<trainData.rows;i++)
			{
				if(trainData.dataTag[i][0]==2)
				{
					dataPositive=dataPositive+enc_traindata[i];
					positiveCount++;
				}
				else
				{
					dataMinus=dataMinus+enc_traindata[i];
					minusCount++;
				}
			}
			x0 = minusCount*dataPositive + positiveCount*dataMinus;
			y0 = minusCount*dataPositive - positiveCount*dataMinus;
			dataPositive = dataPositive * minusCount;
			dataMinus = dataMinus * positiveCount;
			theW = ZZ(2)*minusCount*positiveCount*(dataPositive - dataMinus);

		}

		void testLMC(data testData)
		{
			//encrypt test data
			vec_ZZ encryptedData[100];
			for(int i=0;i<testData.rows;i++)
			{
				encryptedData[i]=encrypt(T, testData.dataSet[i]);
			}

			mat_ZZ M0;
			M0 = innerProdClient(T,S);

			vec_ZZ c0;
			vec_ZZ c;
			ZZ c_star;


			c0 = innerProd(y0,x0,M0);//encrypted data innerMutiply

			vec_ZZ resultC_star;
			resultC_star = decrypt(S,c0);
			cout<<"c_star:"<<resultC_star[0]<<endl;


			c_star=c0[0];
			vec_ZZ x;
			x.SetLength(31);

			int flag=0;

			vec_ZZ enc_tag;
			enc_tag.SetLength(100);

			vec_ZZ temp_enc_tag;
			temp_enc_tag.SetLength(31);

			vec_ZZ temp_dec_tag;
			temp_dec_tag.SetLength(30);

			for(int i=0;i<testData.rows;i++)
			{

				x=encryptedData[i];
				f = theW*x - c_star;
				enc_tag[i] = f; 
				
			}
			//cout<<enc_tag<<endl;
			int count=0;
			int tmp=0;
			int errCount=0;
			while(count!=90)
			{
				for(int i=0;i<31;i++)
				{
					if(i<30)
					{

						temp_enc_tag[i]=enc_tag[count];
						count++;
					}
					else
					{
						temp_enc_tag[i]=ZZ(0);
					}
					
				}
				//cout<<temp_enc_tag<<endl;
				temp_dec_tag=decrypt(S,temp_enc_tag);
				//cout<<temp_dec_tag<<endl;
				for(int j=0;j<30;j++)
				{
					if(temp_dec_tag[j]<ZZ(0))
					{
						flag=4;
					}
					else
					{
						flag=2;
					}
					tmp=count-30+j;
					if(flag!=testData.dataTag[tmp][0])
					{
						errCount++;
					}
					//cout<<flag<<":"<<testData.dataTag[tmp][0]<<endl;
				}	
			}

			float errRate;
			errRate=float(errCount)/float(testData.rows);
			cout<<"ciphertext errRate: "<<errRate<<endl;
		}
		


};



int main() 
{
	data trainData0;
	trainData0 = loadDataSet("testSet.txt",90);
	data testData0;
	testData0 = loadDataSet("testSet.txt",100);


	data trainData;
	data testData;

	trainData.dataSet.SetDims(469,30);
	trainData.dataTag.SetDims(469,1);
	trainData.rows=469;

	testData.dataSet.SetDims(100,30);
	testData.dataTag.SetDims(100,1);
	testData.rows=100;



	mat_RR traindata;//训练数据集 
	traindata.SetDims(469,30);
	mat_RR testdata;
	testdata.SetDims(100,30);


	vec_RR trainTag;//读入训练数据集的标签 
	trainTag.SetLength(469);
	vec_RR testTag;//读入训练数据集的标签 
	testTag.SetLength(100);


	freopen("D_train.txt", "r", stdin);
	for(int i=0; i<469;i++)
	{
		for(int j=0; j<30; j++)
		{
			cin>>traindata[i][j];
			trainData.dataSet[i][j]=TruncToZZ(traindata[i][j]*pow(10,10));
			
		}
	}

	freopen("D_train_result.txt", "r", stdin);
	for(int i=0; i<469;i++)
	{
		
		cin>>trainTag[i];
		trainData.dataTag[i][0]=TruncToZZ(trainTag[i]);
		
	}

	freopen("D_test.txt", "r", stdin);
	for(int i=0; i<100;i++)
	{
		for(int j=0; j<30; j++)
		{
			cin>>testdata[i][j];
			testData.dataSet[i][j]=TruncToZZ(testdata[i][j]*pow(10,10));
		}

	}



	fstream fp_S("D_test_result.txt",fstream::in) ;
	assert(fp_S.is_open());
	for(int i=0;i<100;i++)
	{	
		fp_S >> testTag[i];	
		testData.dataTag[i][0]=TruncToZZ(testTag[i]);
	}
	fp_S.close();


	lmcOnClearText DFLMC(trainData);
  	DFLMC.testLMC(testData);
	lmcOnCiphertext2 CDFLMC2(trainData);
	CDFLMC2.testLMC(testData);
  	




}

