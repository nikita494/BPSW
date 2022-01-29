#include <stdio.h>
#include <gmp.h>


bool miller_rabin(mpz_t n_in, mpz_t b_in){
  mpz_t n, b;
  mpz_init_set(n, n_in);
  mpz_init_set(b, b_in);
  if(mpz_cmp_ui(n, 2) == 0){
    return true;
  }
  else if (mpz_cmp_ui(n, 2) < 0) || mpz_even_p(n)){
    return false;
  }

  if(mpz_cmp_ui(b, 2) < 0){
    b = 2;
  }

  mpz_t g;
  mpz_init(g);
  mpz_gcd(g, n, b);
  while(mpz_cmp_ui(g, 1) != 0){
    if(mpz_cmp(n, g) > 0){
      return false;
    }
    mpz_add_ui(b, b, 1);
    mpz_gcd(g, n, b);
  }

  mpz_t n_1;
  mpz_set(n_1, n);
  mpz_sub_ui(n_1, n_1, 1);

  mpz_t p, q;
  mpz_init(p);
  while(mpz_even_p(n_1)){
    mpz_add_ui(p, p, 1);
    mpz_divexact_ui(n_1, n_1, 2);
  }
  mpz_init_set(q, n_1);

  mpz_t n_1;
  mpz_set(n_1, n);
  mpz_sub_ui(n_1, n_1, 1);

  mpz_t rem;
  mpz_init(rem);
  mpz_powm(rem, b, q, n);
  if(mpz_cmp_ui(rem, 1) == 0 || mpz_cmp(rem, n_1) == 0){
    return true;
  }

  mpz_t i;
  mpz_init_set_ui(i, 1);
  for(mpz_init_set_ui(i, 1); mpz_cmp(i, p) < 0; mpz_add_ui(i, 1)){
    mpz_powm_ui(rem, rem, 2, n);
    if(mpz_cmp(rem, n_1) == 0){
      return true;
    }
  }
  return false;
}
