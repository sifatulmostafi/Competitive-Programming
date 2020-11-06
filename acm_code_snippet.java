/******************************************
 * 
 *	Type: Class File
 * 	Programming Language: Java
 * 	
 * 	Algorithm Used:
 * 	--> Sieve of Atkins
 * 	--> Sieve Of Eratosthenes
 * 	--> Naive Primarility
 * 	--> Fermat Primarility
 * 	--> Millar Rabin Primarility
 *      --> Solovay Strassen Primarility
 *      --> Brute Force Pattern Match
 *      --> Knuth Morris Pratt Pattern Match
 *      --> Floyd Warshall Path Finder
 * 	
 * 	Functions:
 * 	--> GCD
 * 	--> LCM
 * 	--> amountOfDivisors
 * 	--> amountOfFactors
 *  	--> isPrime_naive
 *   	--> isPrime_fermat
 *    	--> isPrime_millarRabin
 * 	--> isPrime_solovayStrasses
 * 	--> jacobi
 * 	--> modPow
 *      --> pow
 *      --> sieve_of_eratosthenes
 *      --> sieve_of_atkins
 *      --> amountOfPrimes
 * 	--> LPD
 * 	--> nCr
 * 	--> largestDigitFinder
 * 	--> isPerfectSquare
 * 	--> baseConverter
 * 	--> lastDigitOfPow
 * 	--> multiplier
 *      --> bin2dec
 *      --> dec2bin
 * 	--> bin2dec_BigIngteger
 * 	--> dec2bin_BigInteger
 *      --> bigSqrt
 * 
 * 	Programming Team: BRACU_Exceptions
 * 	Contributors: Sifatul Mostafi
 * 
 ******************************************/

import java.util.Random;
import java.math.BigInteger;

public class method {

	/*********************************************/
	/** METHOD: gcd (int a,int b) **/
	// returns GCD (Greatest Common Divisor)
	// returns HCF (Highest Common Factor)
	// returns GCF (Highest Common Factor)
	// returns int value
	public static int gcd(int a, int b) {
		int temp;
		while (a > 0) {
			temp = a;
			a = b % a;
			b = temp;
		}
		return b;
	}

	/***********************************************/
	/** METHOD: gcd (long x,long y) **/
	// returns GCD (Greatest Common Divisor)
	// returns HCF (Highest Common Factor)
	// returns GCF (Highest Common Factor)
	// an alternative method of gcd(int a,int b)
	// takes two long values through parameter
	// returns long value
	public static long gcd(long x, long y) {
		if (x == 0)
			return y;
		while (y != 0) {
			if (x > y)
				x = x - y;
			else
				y = y - x;
		}
		return x;
	}

	/***********************************************/
	/** METHOD: GCD(int a, int b) **/
	// returns GCD (Greatest Common Divisor)
	// returns HCF (Highest Common Factor)
	// returns GCF (Highest Common Factor)
	// recursive fuction
	// returns int value
	public static int GCD(int a, int b) {
		return (a == 0) ? b : GCD(b % a, a);
	}

	/***********************************************/
	/** METHOD: LCM(int a,int b) **/
	// returns LCM (Least Common Multiple)
	// uses METHOD GCD(int a,int b) to produce LCD
	// returns int value
	public static int LCM(int a, int b) {
		return (a * b) / GCD(a, b);
	}

	/***********************************************/
	/** METHOD: amountOfdivisors(long num) **/
	// find out the amount of divisors of a number
	// returns long value
	public static int amountOfDivisors(int num) {
		int amountOfDivisor = 0;
		long rootOfNum = (long) Math.sqrt(num) + 1;
		for (int i = 1; i < rootOfNum; i++) {
			if (num % i == 0) {
				amountOfDivisor++;
				if (num / i != i) {
					// System.out.print(num/i+" ");
					amountOfDivisor++;
				}
			}
		}
		return amountOfDivisor;
	}

	/********************************************/
	/** METHOD: amountOfFactors(int n) **/
	// finds the "number of divisors must be needed" to form the number n
	// returns amount of factors
	// returns integer value
	public static int amountOfactors(int n) {
		int factors = 0;
		for (int i = 2; i * i <= n; i++) {
			if (n % i == 0) {
				while (n % i == 0 && n / i > 0) {
					n /= i;
					factors++;
				}
			}
		}
		if (n != 1)
			factors++;
		return factors;
	}

	/**********************************************/
	/** METHOD: isPrime_naive(int num) **/
	// implements naive algorithm
	// determines wether a number prime or not
	// returns boolean value
	public static boolean isPrime_naive(int num) {
		if (num == 2)
			return true;
		if (num % 2 == 0)
			return false;
		int rootOfNum = (int) Math.sqrt(num) + 1;
		for (int i = 2; i < rootOfNum; i += 2)
			if (num % i == 0)
				return false;
		return true;
	}

	/**********************************************/
	/** METHOD: isPrime_fermat(int num) **/
	// implements Fermat primarility algorithm
	// determines wether a number prime or not
	// returns boolean value
	public static boolean isPrime_fermat(long n) {
		// base case
		if (n == 0 || n == 1)
			return false;
		// base case - 2 is prime
		if (n == 2)
			return true;
		// an even number other than 2 is composite
		if (n % 2 == 0)
			return false;
		Random rand = new Random();
		for (int i = 0; i < 50; i++) {
			long r = Math.abs(rand.nextLong());
			long a = r % (n - 1) + 1;
			if (modPow(a, n - 1, n) != 1)
				return false;
		}
		return true;
	}

	/******************************************************/
	/** METHOD: isPrime_millarRabin(int n,int iteration) **/
	// implements Millar Rabin primarility algorithm
	// determines wether a number prime or not
	// returns boolean value
	public static boolean isPrime_millarRabin(long n, int iteration) {
		/** base case **/
		if (n == 0 || n == 1)
			return false;
		/** base case - 2 is prime **/
		if (n == 2)
			return true;
		/** an even number other than 2 is composite **/
		if (n % 2 == 0)
			return false;
		long s = n - 1;
		while (s % 2 == 0)
			s /= 2;
		Random rand = new Random();
		for (int i = 0; i < iteration; i++) {
			long r = Math.abs(rand.nextLong());
			long a = r % (n - 1) + 1, temp = s;
			long mod = modPow(a, temp, n);
			while (temp != n - 1 && mod != 1 && mod != n - 1) {
				// mod = mulMod(mod, mod, n);
				mod = BigInteger.valueOf(mod).pow(2).mod(BigInteger.valueOf(n))
						.longValue();
				temp *= 2;
			}
			if (mod != n - 1 && temp % 2 == 0)
				return false;
		}
		return true;
	}

	/**********************************************/
	/** Function to check if prime or not **/
	// implements Solvay Strassen primarility algorithm
	// determines wether a number prime or not
	// returns boolean value
	public boolean isPrime_SolovayStrassen(long n, int iteration) {
		/** base case **/
		if (n == 0 || n == 1)
			return false;
		/** base case - 2 is prime **/
		if (n == 2)
			return true;
		/** an even number other than 2 is composite **/
		if (n % 2 == 0)
			return false;
		Random rand = new Random();
		for (int i = 0; i < iteration; i++) {
			long r = Math.abs(rand.nextLong());
			long a = r % (n - 1) + 1;
			long jacobian = (n + Jacobi(a, n)) % n;
			long mod = modPow(a, (n - 1) / 2, n);
			if (jacobian == 0 || mod != jacobian)
				return false;
		}
		return true;
	}

	/************************************************/
	/** METHOD: jacobi(long a,long b) **/
	// Function to calculate jacobi (a/b)
	// returns long value
	public static long Jacobi(long a, long b) {
		if (b <= 0 || b % 2 == 0)
			return 0;
		long j = 1L;
		if (a < 0) {
			a = -a;
			if (b % 4 == 3)
				j = -j;
		}
		while (a != 0) {
			while (a % 2 == 0) {
				a /= 2;
				if (b % 8 == 3 || b % 8 == 5)
					j = -j;
			}
			long temp = a;
			a = b;
			b = temp;
			if (a % 4 == 3 && b % 4 == 3)
				j = -j;
			a %= b;
		}
		if (b == 1)
			return j;
		return 0;
	}

	/***********************************************/
	/** METHOD: modPow(long a, long b, long c) **/
	// calculates (a ^ b) % c
	// returns long value
	public static long modPow(long a, long b, long c) {
		long res = 1;
		for (int i = 0; i < b; i++) {
			res *= a;
			res %= c;
		}
		return res % c;
	}

	/**********************************************/
	/** METHOD: pow (long num,long exp) **/
	// returns the power of a Integer/long number upto power exp
	// returns long value
	public static long pow(long num, long exp) {
		long pow = 1;
		while (exp-- > 0)
			pow *= num;
		return pow;
	}

	/***********************************************/
	/** METHOD: sieve_of_eratosthenes(int N) **/
	// returns a boolean array of size of N
	// the array's prime indexed values are only true
	// implements Sieve of Eratosthenes
	public static boolean[] sieve_of_eratosthenes(int N) {
		// initially assume all integers are prime
		boolean[] isPrime = new boolean[N + 1];
		for (int i = 2; i <= N; i++)
			isPrime[i] = true;
		// mark non-primes <= N using Sieve of Eratosthenes
		for (int i = 2; i * i <= N; i++) {
			// if i is prime, then mark multiples of i as nonprime
			// suffices to consider mutiples i, i+1, ..., N/i
			if (isPrime[i]) {
				for (int j = i; i * j <= N; j++) {
					isPrime[i * j] = false;
				}
			}
		}
		return isPrime;
	}

	/***********************************************/
	/** METHOD: sieveOfAtkins(int limit) **/
	// produce primes within the limit
	// prints all the primes within the limit
	public static void sieve_of_atkins(int limit) {
		// initialize the sieve
		boolean[] prime = new boolean[limit + 1];
		prime[2] = true;
		prime[3] = true;
		int root = (int) Math.ceil(Math.sqrt(limit));
		// put in candidate primes: integers which have an odd number of
		// representations by certain quadratic forms
		for (int x = 1; x < root; x++) {
			for (int y = 1; y < root; y++) {
				int n = 4 * x * x + y * y;
				if (n <= limit && (n % 12 == 1 || n % 12 == 5))
					prime[n] = !prime[n];
				n = 3 * x * x + y * y;
				if (n <= limit && n % 12 == 7)
					prime[n] = !prime[n];
				n = 3 * x * x - y * y;
				if ((x > y) && (n <= limit) && (n % 12 == 11))
					prime[n] = !prime[n];
			}
		}
		// eliminate composites by sieving, omit multiples of its square
		for (int i = 5; i <= root; i++) {
			if (prime[i]) {
				for (int j = i * i; j < limit; j += i * i) {
					prime[j] = false;
				}
			}
		}
		System.out.print("\nPrimes = ");
		for (int i = 2; i < prime.length; i++) {
			if (prime[i]) {
				System.out.print(i + " ");
			}
		}
		System.out.println();
	}

	/**********************************************/
	/** METHOD: amountOfPrimes(int N) **/
	// determines and returns the amount of prime numbers within a range
	// uses sieve_of_eratosthenes(int N) method
	// returns the prime amount in integer
	public static int amountOfPrimes(int N) {
		// iplementing sieve_of_eratosthenes
		boolean isPrime[] = sieve_of_eratosthenes(N);
		int primes = 0;
		for (int i = 2; i <= N; i++)
			if (isPrime[i])
				primes++;
		return primes;
	}

	/***********************************************/
	/** METHOD: LPD(long sum ) **/
	// returns LPD (Largest Prime Divisor)
	// returns LPF (Largest Prime Factor)
	public static long LPD(long num) {
		// divide out all of the 2 prime factors
		// use bitwise operators to test even and divide num by 2 for speed
		while ((num & 1L) == 0)
			num = num >> 1;
		// now only working with prime factors > 2 (odd numbers)
		// 3 is the next largest prime factor after 2
		long divisor_current = 3L;
		long lpd = 2L;
		while (divisor_current <= num) {
			if (num % divisor_current == 0L) {
				lpd = divisor_current;
				num = num / divisor_current;
				divisor_current = 3L;
			}
			divisor_current += 1L;
		}
		return lpd;
	}

	/***********************************************/
	/** METHOD: nCr(long n,long r) **/
	// determines the nCr value for n and r
	// returns long value
	public static long nCr(long n, long r) {
		long res = 1, i;
		if (n - r < r)
			r = n - r;
		for (i = 1; i <= r; i++, n--) {
			res = res * n;
			res = res / i;
		}
		return res;
	}

	/***********************************************/
	/** METHOD: largestDigitFinder(String number) **/
	// returns the largest digit in a string
	public static int largestDigitFinder(String number) {
		int largestDigit = 0;
		for (int counter = 0, digit = 0; counter < number.length(); counter++) {
			digit = number.charAt(counter) - 48;
			if (digit > largestDigit)
				largestDigit = digit;
			if (largestDigit == 9)
				return 9;
		}
		return largestDigit;
	}

	/***********************************************/
	/** METHOD: isPerfectSquare(int number) **/
	// determine wether a number is a perfect square or not
	// retuns boolean value
	public static boolean isPerfectSquare(int number) {
		return (Math.sqrt(number) % 1 == 0) ? true : false;
	}

	/***********************************************/
	/** METHOD: isPerfectSquare(long n) **/
	// alternative for the isPerfectSquare(int n)
	// takes long value through parameter
	// returns long value
	public static long isPerfectSquare(long n) {
		int sq = (int) Math.pow(n, .5);
		return (sq * sq == n) ? sq : n;
	}

	/**********************************************/
	/** METHOD: lastDigitOfPow(int base,int exp) **/
	// returns the last digit of base^exp
	// returns int value
	static int lastDigitOfPow(int base, int exp) {
		if (base == 0)
			return 0; // 0^whatever = 0;
		if (exp == 0)
			return 1; // whatever^0 = 1;
		int lastDigit = base;
		while (exp-- > 1) {
			lastDigit *= base;
			lastDigit %= 10;
		}
		return lastDigit % 10;
	}

	/*****************************************************/
	/** METHOD: multiplier(String s1, String s2) **/
	// takes two Integer numbers in String format
	// returns the multiplied value of that two numbers
	public static String multiplier(String s1, String s2) {
		int result[] = new int[s1.length() + s2.length()];
		// iterating over first String
		for (int i = s1.length() - 1; i > -1; i--) {
			// iterating over second String
			for (int j = s2.length() - 1; j > -1; j--) {
				// multiplying value and saving them in result array
				result[i + j + 1] += (s1.charAt(i) - 48) * (s2.charAt(j) - 48);
			}
		}
		// processing the carry operation
		for (int i = result.length - 1; i > -1; i--) {
			if (result[i] > 9) {
				result[i - 1] += result[i] / 10;
				result[i] %= 10;
			}
		}
		// generating a String from the result array
		boolean checker = true;
		String s = "";
		for (int i = 0; i < result.length; i++) {
			if (result[i] != 0 && checker)
				checker = false;
			if (!checker)
				s += result[i];
		}
		if (s.length() == 0)
			return "0";
		return s;
	}

	/************************************************************/
	/** METHOD: baseConverter(String number,int from, int to) **/
	// converts a number from one base to another base
	// returns Stirng value
	public static String baseConverter(String number, int from, int to) {
		return new StringBuilder(Integer.toString(
				Integer.parseInt(number, from), to)).toString();
	}

	/******************************************************************/
	/** METHOD: longDecimalToBase_N (long decimalValue, int Base_N) **/
	// converts a long integer into another base
	// returns Stirng value
	public String longDecimalToBase_N(long decimalValue, int Base_N) {
		return Long.toString(decimalValue, Base_N);
	}

	/******************************************************************/
	/** METHOD: longDecimalToBase_N(String decimalValue,int Base_N) **/
	// converts a string of long integer into another base
	// returns long value
	public long longDecimalToBase_N(String decimalValue, int Base_N) {
		return Long.parseLong(decimalValue, Base_N);
	}

	/******************************************/
	/** METHOD: dec2bin_String (String dec) **/
	// converts a string of long integer into binary
	// returns Stirng value
	public static String dec2bin(String dec) {
		return Long.toBinaryString(Long.parseLong(dec));
	}

	/**********************************/
	/** METHOD: bin2dec(String bin) **/
	// converts a binary number to long
	// returns Stirng value
	public static String bin2dec(String bin) {
		return Long.valueOf(Long.parseLong(bin, 2)).toString();
	}

	/*********************************************/
	/** METHOD: dec2bin_BigInteger(String dec) **/
	// converts a string of decimal value to binary
	// uses BigInteger class
	// returns stirng value
	public static String dec2bin_BigInteger(String dec) {
		return new BigInteger(dec).toString(2);
	}

	/**********************************************/
	/** METHOD: bin2dec_BigInteger (String bin) **/
	// converts a string of binary value to decimal value
	// uses BigInteger class
	// returns string value
	public static String bin2dec_BigInteger(String bin) {
		return new BigInteger(bin, 2).toString();
	}

	/**************************************/
	/** METHOD: bigSqrt(BigInteger num) **/
	// finds out the sqrt() for BigInteger num
	// returns BigInteger value
	static BigInteger bigSqrt(BigInteger num) {
		BigInteger y = num;
		BigInteger div = new BigInteger("2");
		BigInteger x;
		do {
			x = y;
			y = x.add(num.divide(x)).divide(div);
		} while (x.compareTo(y) != 0);
		return y;
	}

	/******************************************************/
	/** METHOD: bigSqrt_fast(BigInteger A) **/
	// another faster implemnatation for BigInteger Sqrt
	// returns BigInteger value
	static BigInteger bigSqrt_fast(BigInteger num) {
		BigInteger temp = num.shiftRight(BigInteger.valueOf(num.bitLength())
				.shiftRight(1).intValue());
		BigInteger result = null;
		while (true) {
			result = temp.add(num.divide(temp)).shiftRight(1);
			if (!temp.equals(result))
				temp = result;
			else
				break;
		}
		return result;
	}

	/***************************************************************/
	/** METHOD: patternPosition_bruteForce(String pat,String txt) **/
	// return offset of first match or N if no match
	// Reads in two strings --> the pattern and the input text
	// searches for the pattern in the input text using brute force.
	public static int patternPosition_bruteForce(String pat, String txt) {
		int M = pat.length();
		int N = txt.length();
		for (int i = 0; i <= N - M; i++) {
			int j;
			for (j = 0; j < M; j++) {
				if (txt.charAt(i + j) != pat.charAt(j))
					break;
			}
			if (j == M)
				return i; // found at offset i
		}
		return N; // not found
	}

	/***************************************************************/
	/** METHOD: patternPosition_bruteForce(char[] pattern,char[] text) **/
	// return offset of first match or N if no match
	// Reads in two strings --> the pattern and the input text
	// searches for the pattern in the input text using brute force.
	public static int search2(char[] pattern, char[] text) {
		int M = pattern.length;
		int N = text.length;
		int i, j;
		for (i = 0, j = 0; i < N && j < M; i++) {
			if (text[i] == pattern[j])
				j++;
			else {
				i -= j;
				j = 0;
			}
		}
		if (j == M)
			return i - M; // found
		else
			return N; // not found
	}

	/*********************************************************************/
	/** METHOD: knuth_morris_pratt() **/
	// return offset of first match or N if no match
	// Reads in two strings --> the pattern and the input text
	// searches for the pattern in the input text using brute force.
	public static int knuth_morris_pratt(String target, String pattern) {
		if (pattern.length() < 2)
			return target.indexOf(pattern.charAt(0)) > 0 ? 1 : 0;
		int[] overlay = new int[pattern.length()];
		overlay[0] = -1;
		overlay[1] = 0;
		int i = 0, j = 1;
		while (j + 1 < pattern.length()) {
			if (pattern.charAt(i) == pattern.charAt(j)) {
				if (i == 0) {
					overlay[j + 1] = 1;
				} else {
					overlay[j + 1] = overlay[j] + 1;
				}
				i++;
				j++;
			} else if (pattern.charAt(j) == pattern.charAt(0)) {
				i = 0;
			} else {
				j++;
			}
		}
		int l = 0, count = 0;
		for (int k = 0; k < target.length(); k++) {
			if (target.charAt(k) == pattern.charAt(l)) {
				if (l == pattern.length() - 1) {
					l = 0;
					count++;
				} else {
					l++;
				}
			} else {
				l = overlay[l] == -1 ? 0 : overlay[l];
			}
		}
		return count;
	}

	/**********************************************************/
	/** METHOD: floyd_warshall(int[] x,int[] y) **/
	// finds minimum necessary steps require over all possible paths between two points.
	// takes two arrays of x co-ordinates and y co-ordinates of all possible paths
	public static double floyd_warshall(int[] x, int[] y) {
		int stones = x.length;
		double d[][] = new double[stones][stones];
		for (int i = 0; i < stones; i++) {
			for (int j = i + 1; j < stones; j++) {
				d[i][j] = d[j][i] = Math.sqrt(Math.pow((x[i] - x[j]), 2)
						+ Math.pow((y[i] - y[j]), 2));
			}
		}
		for (int i = 0; i < stones; i++) {
			for (int j = 0; j < stones; j++) {
				for (int k = 0; k < stones; k++) {
					d[j][k] = Math.min(d[j][k], Math.max(d[j][i], d[i][k]));
				}
			}
		}
		return d[0][1];
	}

	/***********************************************/
	/** METHOD: main(String[] args) **/
	// main method to check other method's functionality
	public static void main(String args[]) {
		// your code to check methods
	}

}
