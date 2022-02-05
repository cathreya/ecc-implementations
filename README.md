This repo has my implementations of Reed-Solomon codes (and their efficient decoding and list-decoding) as well as Reed-Muller codes (and their efficient decoding). 

## Reed-Solomon Codes
The description of the various requirements can be found "Programming_Assignment_1__RS_Codes_and_List_Decoding_of_RS_Codes.pdf". The code implements efficient decoding using the Berlekamp-Welch Algorithm. It also implements List decoding based on similar principals.

### Dependencies
- numpy 1.20.1 + pandas 1.2.3 - for visualizing the matrix during Gaussian Elimination.
- sympy 1.7.1 - used for manipulating polynomials in F_2, like addition, multiplication, inverse etc. Sympy [only supports polynomials over prime fields](https://github.com/sympy/sympy/issues/21104). There are some hacks that allow some operations (for example- declaring a finite extension of F_2), however these are also extremely limited. Therefore to perform operations over F_{2^m} as was required in the assignment, I had to write my own functions (such as for division). However I was unable to implement multivariate factorization.
- sage 9.2 - only for multivariate factorization used in Q4/Q5. Sage is only used in a separate script called "factorize_using_sage". I am not entirely sure whether the function is reliable since I couldn't verify by hand in situations where it failed.


### Status of work

- Q1,2,3 have been implemented completely. As far as my testing goes, these seem to work correctly.

- Q4, the interpolation stage is implemented correctly and verified by evaluating Q(x,y) at alpha_i, y_i for each i. Also it satisfies the degree constraint by definition and is non zero. However the factorization stage does not seem to work as expected (even when there are no errors introduced in the codeword). For example the input:

```
Enter n 10
Enter k 4
Enter the number of non-zero coefficients 2
Enter power of X 0
Enter p, where coefficient in the form beta^p (beta is the primitive element) 0
Enter power of X 1
Enter p, where coefficient in the form beta^p (beta is the primitive element) 2
```

works as intended and produces the correct output. However, 

```
Enter n 10
Enter k 4
Enter the number of non-zero coefficients 3
Enter power of X 0
Enter p, where coefficient in the form beta^p (beta is the primitive element) 1
Enter power of X 1
Enter p, where coefficient in the form beta^p (beta is the primitive element) 2
Enter power of X 2
Enter p, where coefficient in the form beta^p (beta is the primitive element) 3
```

returns a list containing a single all zero codeword (which is incorrect).

- Q5, the interpolation stage is implemented correctly and verified by evaluating Q(x,y) at alpha_i, y_i for each i. It also functions correctly when there are no errors. 
With a single error it usually returns a correct list. However, if there is more than one error, the polynomial factorizaton doesn't go through and an empty list is returned.


### Learnings
- My foremost learning in this assignment has been about working with non-prime finite fields. Previous to this, I had only worked with prime finite fields. The assignment helped me understand the parallel between manipulating polynomials modulo an irreducible polynomial and integers modulo a prime. 
- The assignment has also helped me understand the gap between theory and practice. I would assume to be useful in real life, these decoding algorithms must be run several times, and so must be implemented much more efficiently than I have here. For example: operations like matrix vector multiplication is heavily optimized for integers. However, since Sympy doesn't support non-prime fields, I had to write my own function which was significantly slower.
- I've also learnt not to blindly trust (even reputable) libraries. Several function calls return an incorrect output rather than a NotImplemented exception. In the course of the assignment, I have discovered and reported 3 bugs. [2 in sympy](https://github.com/sympy/sympy/issues/21104): a) sympy reported a field element to not have an exact inverse and b) The factor() function over Finite Extensions simply returns the function as is. [1 in sage](https://ask.sagemath.org/question/56220/finite-field-extension/): a programming bug when using finite extensions. 

## Reed-Muller codes
The description of the various requirements can be found "Programming_Assignment_2__RM_Codes_and_Decoding.pdf". The code implements efficient decoding using the Majority Logic decoding algorithm. It also implements a decoding algorithm based on the Fast Hadamard transform when the parameters are appropriate.

### Dependencies
- numpy 1.20.1 - for pretty printing arrays and verifying operations
- sympy 1.7.1 - used for manipulating polynomials in F_2

### Status of work
- All problems have been implemented completely. However Q1(d) takes a lot of time to execute since it requires O(2^20 (2^20 positions) * 20 (variables) * 10 (monomials)) = 209715200 operations. Similarly, for performing MLD decoding we require O(2^20 * 616666) which is 646621167616 operations. Hence this Q2(d) does not terminate.

### Learnings
- Both decoding algorithms can be succinctly implemented. In fact most of the implementation didn't even require the use of any library and can be coded up even in lower level languages like C++ quite easily. This made me appreciate these RM codes from a practical point of view.