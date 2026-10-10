#include <bits/stdc++.h>
using namespace std;

// ----- term -----
// Represents a * x^b * sin^c(x) * cos^d(x)
struct term {
    long long a; // coefficient
    int b, c, d; // exponents for x, sin, cos
    term() : a(0), b(0), c(0), d(0) {}
    term(long long _a, int _b, int _c, int _d) : a(_a), b(_b), c(_c), d(_d) {}
    bool same_key(const term &o) const { return b==o.b && c==o.c && d==o.d; }
};

// ----- poly -----
struct poly {
    vector<term> t; // list of terms
    poly() {}
    explicit poly(const term &single) { if(single.a!=0) t.push_back(single); }
    // combine like terms and remove zero coefficients
    void simplify() {
        sort(t.begin(), t.end(), [](const term &x, const term &y){
            if (x.b!=y.b) return x.b>y.b; // descending x exponent
            if (x.c!=y.c) return x.c>y.c; // then sin exponent
            return x.d>y.d; // then cos exponent
        });
        vector<term> res;
        for (auto &tr: t) {
            if (tr.a==0) continue;
            if (!res.empty() && res.back().same_key(tr)) {
                res.back().a += tr.a;
            } else {
                res.push_back(tr);
            }
        }
        // remove zero after merging
        t.clear();
        for (auto &tr: res) if (tr.a!=0) t.push_back(tr);
    }
    // addition
    poly operator+(const poly &o) const {
        poly ans;
        ans.t.reserve(t.size()+o.t.size());
        ans.t.insert(ans.t.end(), t.begin(), t.end());
        ans.t.insert(ans.t.end(), o.t.begin(), o.t.end());
        ans.simplify();
        return ans;
    }
    // subtraction
    poly operator-(const poly &o) const {
        poly ans;
        ans.t.reserve(t.size()+o.t.size());
        ans.t.insert(ans.t.end(), t.begin(), t.end());
        for (auto &tr: o.t) {
            term neg = tr; neg.a = -neg.a; ans.t.push_back(neg);
        }
        ans.simplify();
        return ans;
    }
    // multiplication
    poly operator*(const poly &o) const {
        poly ans;
        if (t.empty() || o.t.empty()) return ans; // zero
        ans.t.reserve(t.size()*o.t.size());
        for (auto &l: t) {
            for (auto &r: o.t) {
                term cur;
                cur.a = l.a * r.a;
                cur.b = l.b + r.b;
                cur.c = l.c + r.c;
                cur.d = l.d + r.d;
                ans.t.push_back(cur);
            }
        }
        ans.simplify();
        return ans;
    }
    // derivative of a polynomial (term-wise)
    poly derivate() const {
        poly ans;
        for (auto &tr: t) {
            // derivative of a * x^b * sin^c * cos^d
            // term1: a*b * x^{b-1} * sin^c * cos^d
            if (tr.b != 0) {
                term t1(tr.a * tr.b, tr.b-1, tr.c, tr.d);
                ans.t.push_back(t1);
            }
            // term2: a*c * x^{b} * sin^{c-1} * cos^{d+1}
            if (tr.c != 0) {
                term t2(tr.a * tr.c, tr.b, tr.c-1, tr.d+1);
                ans.t.push_back(t2);
            }
            // term3: -a*d * x^{b} * sin^{c+1} * cos^{d-1}
            if (tr.d != 0) {
                term t3(-tr.a * tr.d, tr.b, tr.c+1, tr.d-1);
                ans.t.push_back(t3);
            }
        }
        ans.simplify();
        return ans;
    }
    // output as string according to problem format
    string to_string_poly() const {
        if (t.empty()) return "0";
        string out;
        for (size_t i=0;i<t.size();++i) {
            const term &tr = t[i];
            long long coeff = tr.a;
            bool negative = coeff < 0;
            long long abscoeff = llabs(coeff);
            // sign handling
            if (i==0) {
                if (negative) out.push_back('-');
            } else {
                out.push_back(negative ? '-' : '+');
            }
            // coefficient output
            bool need_coeff = (abscoeff != 1) || (tr.b==0 && tr.c==0 && tr.d==0);
            if (need_coeff) out += to_string(abscoeff);
            // x part
            if (tr.b>0) {
                out.push_back('x');
                if (tr.b!=1) out += "^" + to_string(tr.b);
            }
            // sin part
            if (tr.c>0) {
                out += "sin";
                if (tr.c!=1) out += "^" + to_string(tr.c);
                out.push_back('x');
            }
            // cos part
            if (tr.d>0) {
                out += "cos";
                if (tr.d!=1) out += "^" + to_string(tr.d);
                out.push_back('x');
            }
        }
        return out;
    }
};

// ----- fraction -----
struct frac {
    poly p, q; // numerator / denominator
    frac() {}
    explicit frac(const poly &num) : p(num), q(poly(term(1,0,0,0))) {}
    explicit frac(const term &t) : p(poly(t)), q(poly(term(1,0,0,0))) {}
    // addition
    frac operator+(const frac &o) const {
        // p/q + o.p/o.q = (p*o.q + o.p*q) / (q*o.q)
        frac ans;
        ans.p = p*o.q + o.p*q;
        ans.q = q*o.q;
        return ans;
    }
    // subtraction
    frac operator-(const frac &o) const {
        frac ans;
        ans.p = p*o.q - o.p*q;
        ans.q = q*o.q;
        return ans;
    }
    // multiplication
    frac operator*(const frac &o) const {
        frac ans;
        ans.p = p*o.p;
        ans.q = q*o.q;
        return ans;
    }
    // division
    frac operator/(const frac &o) const {
        frac ans;
        ans.p = p*o.q;
        ans.q = q*o.p;
        return ans;
    }
    // derivative
    frac derivate() const {
        // (p/q)' = (p'*q - p*q') / (q*q)
        frac ans;
        ans.p = p.derivate()*q - p*q.derivate();
        ans.q = q*q;
        return ans;
    }
    // output
    string to_string_frac() const {
        // simplify zero numerator
        if (p.t.empty()) return "0";
        // denominator 1
        bool denom_one = (q.t.size()==1 && q.t[0].a==1 && q.t[0].b==0 && q.t[0].c==0 && q.t[0].d==0);
        string num = p.to_string_poly();
        if (denom_one) return num;
        string den = q.to_string_poly();
        // need parentheses if multiple terms in numerator or denominator
        if (p.t.size()>1) num = "(" + num + ")";
        if (q.t.size()>1) den = "(" + den + ")";
        return num + "/" + den;
    }
};

// ----- lexer -----
enum TokenType {TOK_END, TOK_NUM, TOK_X, TOK_SIN, TOK_COS, TOK_PLUS, TOK_MINUS, TOK_MUL, TOK_DIV, TOK_POW, TOK_LP, TOK_RP};
struct Token { TokenType type; long long value; };

string s; int pos, n;
Token cur;

void skip_spaces(){ while(pos<n && isspace(s[pos])) ++pos; }
Token next_token(){
    skip_spaces();
    if(pos>=n) return {TOK_END,0};
    char c=s[pos];
    if(isdigit(c)){
        long long val=0; while(pos<n && isdigit(s[pos])){ val=val*10+(s[pos]-'0'); ++pos; }
        return {TOK_NUM,val};
    }
    if(c=='x'){ ++pos; return {TOK_X,0}; }
    if(c=='s' && s.substr(pos,3)=="sin"){ pos+=3; return {TOK_SIN,0}; }
    if(c=='c' && s.substr(pos,3)=="cos"){ pos+=3; return {TOK_COS,0}; }
    ++pos;
    switch(c){
        case '+': return {TOK_PLUS,0};
        case '-': return {TOK_MINUS,0};
        case '*': return {TOK_MUL,0};
        case '/': return {TOK_DIV,0};
        case '^': return {TOK_POW,0};
        case '(': return {TOK_LP,0};
        case ')': return {TOK_RP,0};
        default: return {TOK_END,0};
    }
}

void advance(){ cur = next_token(); }

bool token_starts_factor(TokenType t){
    return t==TOK_NUM || t==TOK_X || t==TOK_SIN || t==TOK_COS || t==TOK_LP;
}

// forward declarations
frac parse_expr();
frac parse_term();
frac parse_factor();

frac parse_expr(){
    frac left = parse_term();
    while(cur.type==TOK_PLUS || cur.type==TOK_MINUS){
        Token op = cur; advance();
        frac right = parse_term();
        if(op.type==TOK_PLUS) left = left + right; else left = left - right;
    }
    return left;
}

// term handles *, / and implicit multiplication
frac parse_term(){
    frac left = parse_factor();
    while(true){
        if(cur.type==TOK_MUL){ advance(); frac right = parse_factor(); left = left * right; }
        else if(cur.type==TOK_DIV){ advance(); frac right = parse_factor(); left = left / right; }
        else if(token_starts_factor(cur.type)){
            // implicit multiplication
            frac right = parse_factor(); left = left * right;
        } else break;
    }
    return left;
}

frac parse_factor(){
    if(cur.type==TOK_PLUS){ advance(); return parse_factor(); }
    if(cur.type==TOK_MINUS){ advance(); // unary minus
        // represent as -1 * factor
        frac f = parse_factor();
        return frac(term(-1,0,0,0)) * f;
    }
    if(cur.type==TOK_NUM){ long long val = cur.value; advance(); return frac(term(val,0,0,0)); }
    if(cur.type==TOK_X){ advance(); int exp=1; if(cur.type==TOK_POW){ advance(); if(cur.type==TOK_NUM){ exp = (int)cur.value; advance(); } }
        return frac(term(1,exp,0,0)); }
    if(cur.type==TOK_SIN || cur.type==TOK_COS){
        bool is_sin = (cur.type==TOK_SIN);
        advance(); int exp=1; if(cur.type==TOK_POW){ advance(); if(cur.type==TOK_NUM){ exp = (int)cur.value; advance(); } }
        // expect 'x'
        if(cur.type!=TOK_X){ /* error */ }
        advance();
        if(is_sin) return frac(term(1,0,exp,0)); else return frac(term(1,0,0,exp));
    }
    if(cur.type==TOK_LP){ advance(); frac inside = parse_expr(); if(cur.type==TOK_RP) advance(); return inside; }
    // fallback zero
    return frac(term(0,0,0,0));
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    string input; if(!(cin>>input)) return 0;
    s = input; pos=0; n=s.size(); advance();
    frac expr = parse_expr();
    // output original expression as fraction
    cout << expr.to_string_frac() << "\n";
    frac deriv = expr.derivate();
    cout << deriv.to_string_frac() << "\n";
    return 0;
}
