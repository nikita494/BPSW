#include <stdio.h>
#include <stdbool.h>
#include <gmp.h>

void transform_num(mpz_t n_in, mpz_t p, mpz_t q){
    mpz_t n;
    mpz_init_set(n, n_in);
    while(mpz_even_p(n)){
        mpz_add_ui(p, p, 1);
        mpz_divexact_ui(n, n, 2);
    }
    mpz_set(q, n);
    mpz_clear(n);
}

void mulmod(mpz_t res, mpz_t n1, mpz_t n2, mpz_t mod){
    mpz_mul(res, n1, n2);
    mpz_mod(res, res, mod);
}

bool miller_rabin(mpz_t n_in, mpz_t b_in){
    mpz_t n, b;
    mpz_init_set(n, n_in);
    mpz_init_set(b, b_in);
    if(mpz_cmp_ui(n, 2) == 0){
        mpz_clears(n, b, NULL);
        return true;
    }
    else if (mpz_cmp_ui(n, 2) < 0 || mpz_even_p(n)){
        mpz_clears(n, b, NULL);
        return false;
    }
    if(mpz_cmp_ui(b, 2) < 0){
        mpz_set_ui(b, 2);
    }

    mpz_t g;
    mpz_init(g);
    mpz_gcd(g, n, b);
    while(mpz_cmp_ui(g, 1) != 0){
        if(mpz_cmp(n, g) > 0){
            mpz_clears(n, b, g, NULL);
            return false;
        }
        mpz_add_ui(b, b, 1);
        mpz_gcd(g, n, b);
    }
    mpz_clear(g);

    mpz_t n_1;
    mpz_init_set(n_1, n);
    mpz_sub_ui(n_1, n_1, 1);

    mpz_t p, q;
    mpz_init(p);
    mpz_init(q);
    transform_num(n_1, p, q);

    mpz_t rem;
    mpz_init(rem);
    mpz_powm(rem, b, q, n);
    if(mpz_cmp_ui(rem, 1) == 0 || mpz_cmp(rem, n_1) == 0){
        mpz_clears(n, b, n_1, p, q, rem, NULL);
        return true;
    }

    mpz_t i;
    for(mpz_init_set_ui(i, 1); mpz_cmp(i, p) < 0; mpz_add_ui(i, i, 1)){
        mpz_powm_ui(rem, rem, 2, n);
        if(mpz_cmp(rem, n_1) == 0){
            mpz_clears(n, b, n_1, p, q, rem, i, NULL);
            return true;
        }
    }
    mpz_clears(n, b, n_1, p, q, rem, i, NULL);
    return false;
}

bool lucas_selfridge(mpz_t n_in){
    mpz_t n;
    mpz_init_set(n, n_in);
    if(mpz_cmp_ui(n, 2) == 0){
        mpz_clear(n);
        return true;
    }
    else if(mpz_cmp_ui(n, 2) < 0 || mpz_even_p(n)){
        mpz_clear(n);
        return false;
    }
    else if(mpz_perfect_square_p(n) == 0){
        mpz_clear(n);
        return false;
    }

    mpz_t gg;
    mpz_init(gg);
    mpz_t abs;
    mpz_init_set_ui(abs, 3);
    mpz_t sign;
    mpz_init_set_ui(sign, -1);
    mpz_t jacobi;
    mpz_init(jacobi);
    do{
        mpz_neg(sign, sign);
        mpz_add_ui(abs, abs, 2);
        mpz_t g;
        mpz_gcd(g, n, abs);
        if(mpz_cmp_ui(g, 1) > 0 & mpz_cmp(g, n) < 0){
            return false;
        }
        mpz_mul(gg, abs, sign);

    }while(mpz_jacobi(gg, n) != 1);
    mpz_clears(abs, sign, jacobi, NULL);

    mpz_t p, q;
    mpz_init_set_ui(p, 1);
    mpz_init_set(q, p);
    mpz_mul(q, p, p);
    mpz_sub(q, gg);
    mpz_divexact_ui(q, q, 4);

    mpz_t n_1, s, d;
    mpz_init_set(n_1, n);
    mpz_add_ui(n_1, n_1, 1);
    mpz_init(s);
    mpz_init(d);
    transform_num(n_1, s, d);

    mpz_t u, v, u2m, v2m, qm, qm2, qkd, bit, bits;
    mpz_init_set_ui(u, 1);
    mpz_init_set(v, p);
    mpz_init_set_ui(u2m, 1);
    mpz_init_set(v2m, p);
    mpz_init_set(qm, q);
    mpz_init_set(qm2, q);
    mpz_mul_ui(qm2, qm2, 2);
    mpz_init_set(qkd, q);
    mpz_init_set_ui(bits, mpz_sizeinbase(d, 2));
    for(mpz_init_set_ui(bit, 1); mpz_cmp(bit, bits) < 0; mpz_add_ui(bit, 1)){
        mulmod(u2m, u2m, v2m, n);
        mulmod(v2m, v2m, v2m, n);
        while(mpz_cmp(v2m, qm2) < 0){
            mpz_add(v2m, v2m, n);
        }
        mpz_sub(v2m, v2m, qm2);
        mulmod(qm, qm, qm, n);
        mpz_set(qm2, qm);
        mpz_mul_ui(qm2, qm2, 2);
        if(mpz_tstbit(d, bit)){
            mpz_t t1, t2;
            mpz_init_set(t1, u2m);
            mulmod(t1, t1, v, n);
            mpz_init_set(t2, v2m);
            mulmod(t2, t2, u, n);

            mpz_t t3, t4;
            mpz_init_set(t3, v2m);
            mulmod(t3, t3, v, n);
            mpz_init_set(t4, u2m);
            mulmod(t4, t4, u, n);
            mulmod(t4, t4, dd, n);

            mpz_add(u, t1, t2);
            if(!mpz_even_p(u)){
                mpz_add(u, u, n);
            }
            mpz_divexact_ui(u, u, 2);
            mpz_mod(u, u, n);

            mpz_add(v, t3, t4);
            if(!mpz_even_p(v)){
                mpz_add(v, v, n);
            }
            mpz_divexact_ui(v, v, 2);
            mpz_mod(v, v, n);
            mulmod(qkd, qkd, qm, n);
        }
    }

    if(mpz_cmp_ui(u, 0) == 0 || mpz_cmp_ui(v, 0) == 0){
        return true;
    }

    mpz_set(qkd2, qkd);
    mpz_mul_ui(qkd2, qkd2, 2);
    mpz_t r;
    mpz_t s_1;
    mpz_init_set(s_1, s);
    mpz_sub_ui(s_1, s_1, 1);
    for(mpz_init_set_ui(r, 1); mpz_cmp(r, s) <  0; mpz_add_ui(r, r, 1)){
        mulmod(v, v, v, n);
        mpz_sub(v, v, qkd2);
        if(mpz_cmp_ui(v, 0) < 0) mpz_add(v, v, n);
        if(mpz_cmp_ui(v, 0) < 0) mpz_add(v, v, n);
        if(mpz_cmp(v, n) >= 0) mpz_sub(v, v, n);
        if(mpz_cmp(v, n) >= 0) mpz_sub(v, v, n);
        if(mpz_cmp_ui(v, 0) == 0){
            return true;
        }
        if(mpz_cmp(r, s_1) < 0){
            mulmod(qkd, qkd, qkd, n);
            mpz_set(qkd2, qkd);
            mpz_mul_ui(qkd2, qkd2, 2);
        }
    }
    return false;
}
