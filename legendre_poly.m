function Pn = legendre_poly(n, x)
% 计算n阶勒让德多项式 Pn(x) 的值
    if n < 0
        error('勒让德多项式的阶数 n 必须为非负整数。');
    end
    Pn = 0;
    k_max = floor(n/2);
    for k = 0:k_max
        numerator = (-1)^k * factorial(2*n - 2*k);
        denominator = 2^n * factorial(k) * factorial(n - k) * factorial(n - 2*k);
        term = numerator / denominator * (x .^ (n - 2*k));
        Pn = Pn + term;
    end
end