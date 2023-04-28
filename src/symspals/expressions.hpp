#ifndef NODEINCLUDED
#define NODEINCLUDED

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include <type_traits>
#include <algorithm>
#include <initializer_list>
#include "common.hpp"
using namespace std;
namespace symspals
{
    class Expression; // forward declaration
    class Sym;        // forward declaration

    class Node
    {
    public:
        virtual ostream &print(ostream &os) { return os << "Node"; };
        virtual bool equals(const Node &node) const { return this == &node; };
        virtual bool is_sym() const { return false; };
        virtual bool is_const() const { return false; };
        virtual bool is_zero() const { return false; };
        virtual bool is_plus() const { return false; };
        virtual bool is_mult() const { return false; };
        virtual bool is_unary() const { return false; };
        virtual bool is_binary() const { return false; };
        virtual bool is_leaf() const { return false; };
        virtual bool is_factor() const = 0;
        virtual bool depends_on(const Expression &sym) const = 0;
        virtual const Expression &dep(const int i) const
        {
            throw runtime_error("Node::dep() not available");
        };
        virtual double value() const
        {
            throw runtime_error("Node::value() not available");
        };
    };

    class Expression : public shared_ptr<Node>
    {
    public:
        Expression() : shared_ptr<Node>(nullptr) {}
        Expression(const double val); // forward declaration, implemented after ConstNode is defined;
        Expression(const shared_ptr<Node> &ptr) : shared_ptr<Node>(ptr){};
        template <typename Derived>
        static Expression make_new(const Derived &node)
        {
            return Expression(make_shared<Derived>(node));
        }
        friend ostream &operator<<(ostream &os, const Expression &expr) { return expr->print(os); };
        // bool operator==(const Expression &expr) const { return this->get()->equals(*expr.get()); };
    };

    class LeafNode : public Node
    {
        bool is_leaf() const override { return true; };
        virtual bool is_factor() const override { return true; };
        virtual bool depends_on(const Expression &sym) const override { return false; };
    };

    class SymNode : public LeafNode
    {
    public:
        SymNode(const string &name) : name(name){};
        const string name;
        bool is_sym() const override { return true; };
        virtual ostream &print(ostream &os) { return os << name; };
        virtual bool depends_on(const Expression &sym) const override; // forward declaration
    };

    class ConstNode : public LeafNode
    {
    public:
        ConstNode(const double value) : value_(value){};
        const double value_;
        bool is_const() const override { return true; };
        bool is_zero() const override { return value_ == 0; };
        double value() const override { return value_; };
        virtual ostream &print(ostream &os) { return os << value_; };
        virtual bool equals(const Node &node) const
        {
            if (!node.is_const())
                return false;
            else
                return value() == node.value();
        };
    };
    Expression::Expression(const double val) : shared_ptr<Node>(make_new<ConstNode>(val)){

                                               };

    class ZeroNode : public LeafNode
    {
    public:
        ZeroNode(){};
        virtual ostream &print(ostream &os) { return os << "00"; };
        virtual bool equals(const Node &node) const
        {
            return node.is_zero();
        };
        bool is_zero() const override { return true; };
        bool is_const() const override { return true; };
        double value() const override { return 0; };
    };

    class UnaryNode : public Node
    {
    public:
        virtual const Expression &dep(const int i) const override { return dep1; };
        bool is_unary() const override { return true; };
        Expression dep1;
        virtual bool is_factor() const override { return dep1->is_factor(); };
        virtual bool depends_on(const Expression &sym) const override; // forward declaration
    };

    class BinaryNode : public Node
    {
    public:
        BinaryNode(const Expression &expr1, const Expression &expr2) : expr1(expr1), expr2(expr2){};
        virtual const Expression &dep(const int i) const override { return (i == 0) ? expr1 : expr2; };
        bool is_binary() const override { return true; };
        const Expression expr1;
        const Expression expr2;
        virtual bool is_factor() const override { return expr1->is_factor() && expr2->is_factor(); };
        virtual bool depends_on(const Expression &sym) const override; // forward declaration
    };

    class PlusNode : public BinaryNode
    {
    public:
        PlusNode(const Expression &expr1, const Expression &expr2) : BinaryNode(expr1, expr2){};
        virtual ostream &print(ostream &os) override { return os << "(" << expr1 << "+" << expr2 << ")"; };
        bool is_plus() const override { return true; };
        virtual bool is_factor() const override { return false; };
    };
    class MultNode : public BinaryNode
    {
    public:
        MultNode(const Expression &expr1, const Expression &expr2) : BinaryNode(expr1, expr2){};
        virtual ostream &print(ostream &os) override { return os << "(" << expr1 << "*" << expr2 << ")"; };
        bool is_mult() const override { return true; };
    };

    class Sym : public Expression
    {
    public:
        Sym(const string &name) : Expression(Expression::make_new(SymNode(name))){};
    };

    bool SymNode::depends_on(const Expression &sym) const
    {
        return this->equals(*sym);
    };
    bool UnaryNode::depends_on(const Expression &sym) const { return dep1->depends_on(sym); };
    bool BinaryNode::depends_on(const Expression &sym) const { return expr1->depends_on(sym) || expr2->depends_on(sym); };

    class Const : public Expression
    {
    public:
        Const(const double value) : Expression(Expression::make_new(ConstNode(value))){};
    };

    // todo make this friend functions of Expression

    Expression operator+(const Expression &expr1, const Expression &expr2)
    {
        // simplifications
        if (expr1->is_zero())
        {
            return expr2;
        }
        if (expr2->is_zero())
        {
            return expr1;
        }
        if (expr1->is_const() && expr2->is_const())
        {
            return expr1->value() + expr2->value();
        }
        return Expression::make_new(PlusNode(expr1, expr2));
    }

    Expression operator*(const Expression &expr1, const Expression &expr2)
    {
        // simplificiations
        if (expr1->is_zero() || expr2->is_zero())
        {
            return 0.0;
        }
        if (expr1->is_const() && expr2->is_const())
        {
            return expr1->value() * expr2->value();
        }
        if (expr1->equals(ConstNode(1.0)))
        {
            return expr2;
        }
        if (expr2->equals(ConstNode(1.0)))
        {
            return expr1;
        }
        return Expression::make_new(MultNode(expr1, expr2));
    }

    Expression operator-(const Expression &expr1)
    {
        return -1.0 * expr1;
    }
    Expression operator-(const Expression &expr1, const Expression &expr2)
    {
        return expr1 + (-1.0 * expr2);
    }

    class ExprVec : public vector<Expression>
    {
    public:
        ExprVec(const vector<Expression> &sym) : vector<Expression>(sym){};
        bool appears_in(const Expression &expr) const
        {
            for (auto &sym : *this)
            {
                if (expr->depends_on(sym))
                {
                    return true;
                }
            }
            return false;
        }
    };

    vector<Expression> GetFactors(const Expression &expr, const ExprVec &syms)
    {
        // check if expr is factorable
        vector<Expression> factors;
        // base case
        if ((expr->is_leaf() || (!syms.appears_in(expr))) || expr->is_unary())
        {
            factors.push_back(expr);
        }
        else if (expr->is_mult())
        {
            if (!((syms.appears_in(expr->dep(0))) ^ (syms.appears_in(expr->dep(1)))))
            {
                throw std::runtime_error("Expression is not factorable - both sides of multiplication contain symbols");
            }
            auto fac1 = GetFactors(expr->dep(0), syms);
            auto fac2 = GetFactors(expr->dep(1), syms);
            // concatenate fac1 and fac2
            factors.insert(factors.end(), fac1.begin(), fac1.end());
            factors.insert(factors.end(), fac2.begin(), fac2.end());
        }
        else
        {
            throw std::runtime_error("Expression is not factorable");
        }
        return factors;
    }

    vector<Expression> GetTerms(const Expression &expr, const ExprVec &syms)
    {
        vector<Expression> terms;
        // check if binary node
        if (expr->is_factor())
        {
            terms.push_back(expr);
        }
        else if (expr->is_binary())
        {
            auto expr1 = expr->dep(0);
            auto expr2 = expr->dep(1);
            // check if plus node
            if (expr->is_plus())
            {
                auto terms1 = GetTerms(expr1, syms);
                auto terms2 = GetTerms(expr2, syms);
                terms.insert(terms.end(), terms1.begin(), terms1.end());
                terms.insert(terms.end(), terms2.begin(), terms2.end());
            }
            else if (expr->is_mult() && (expr1->is_plus() && syms.appears_in(expr1)))
            {
                // (a+b)*c = ac + bc
                auto termss = GetTerms(expr1->dep(0) * expr2 + expr1->dep(1) * expr2, syms);
                terms.insert(terms.end(), termss.begin(), termss.end());
            }
            else if (expr->is_mult() && (expr2->is_plus() && syms.appears_in(expr2)))
            {
                // a*(b+c) = ab + ac
                auto termss = GetTerms(expr1 * expr2->dep(0) + expr1 * expr2->dep(1), syms);
                terms.insert(terms.end(), termss.begin(), termss.end());
            }
            else if (expr->is_mult())
            {
                terms.push_back(expr);
            }
            else
            {
                throw runtime_error("GetTerms: unknown binary node");
            }
        }
        else
        {
            throw runtime_error("GetTerms: unknown node");
        }
        return terms;
    };

    template <typename T>
    class TripletVec : public vector<Triplet<T>>
    {
    public:
        T &get_el(const int i, const int j)
        {
            for (auto &triplet : *this)
            {
                if (triplet.index.row == i && triplet.index.col == j)
                {
                    return triplet.value;
                }
            }
            this->push_back(Triplet<T>(i, j, 0.0));
            return this->back().value;
        }
    };

    TripletVec<Expression> GetCoefficients(const vector<Expression> &expr_vec, const vector<Expression> &sym_vec)
    {
        TripletVec<Expression> ret;
        for (size_t i = 0; i < expr_vec.size(); i++)
        {
            Expression expr = expr_vec.at(i);
            // vector<Expression> coefficients(sym_vec.size(), 0.0);
            // get the terms of the expression
            vector<Expression> terms = GetTerms(expr, sym_vec);
            // iterate over all terms
            for (auto term : terms)
            {
                // iterate over all symbols
                int count = 0;
                int coeff_ind = 0;
                Expression coeff(1.0);
                for (size_t j = 0; j < sym_vec.size(); j++)
                {
                    auto sym = sym_vec.at(j);
                    // check if the symbol is in the term
                    if (term->depends_on(sym))
                    {
                        vector<vector<Expression>> terms_factored;
                        auto factors = GetFactors(term, sym_vec);
                        terms_factored.push_back(factors);
                        coeff_ind = j;
                        for (auto factor : factors)
                        {
                            if (!factor->equals(*sym))
                            {
                                coeff = coeff * factor;
                            }
                        }
                        count++;
                    }
                }
                if (count > 1)
                {
                    throw runtime_error("GetCoefficients: expression is not linear");
                }
                else if (count == 1)
                {
                    ret.get_el(i, coeff_ind) = ret.get_el(i, coeff_ind) + coeff;
                }
            }
        }
        return ret;
    };

    void OrderDepthFirstRecurse(const Expression expr, vector<Expression> &result)
    {
        // // check if expr is already in result
        // if (find(result.begin(), result.end(), expr) != result.end())
        // {
        //     return;
        // }
        // check if expression is a binary
        if (expr->is_binary())
        {
            auto expr1 = expr->dep(0);
            auto expr2 = expr->dep(1);
            OrderDepthFirstRecurse(expr1, result);
            OrderDepthFirstRecurse(expr2, result);
        }
        result.push_back(expr);
    };
    class Function
    {
        // the set of possible instructions is very limited
        enum Instruction
        {
            INPUT,
            OUTPUT,
            PLUS,
            MULT,
            CONST
        };
        struct AlgEl
        {
            Instruction ins;
            int arg1;
            int arg2;
            double val;
        };

    public:
        Function(){};
        Function(const vector<Expression> &input, const vector<Expression> &output)
        {
            // TODO: at this point no re-use of workspace is done, the lifetime of workspace variables is not taken into account
            vector<Expression> ordered_expression;
            for (auto expr : output)
                OrderDepthFirstRecurse(expr, ordered_expression);
            algorithm.reserve(ordered_expression.size());
            work.reserve(ordered_expression.size());
            for (auto expr : ordered_expression)
            {
                // check if expr is a Sym
                if (expr->is_sym())
                {
                    // find the index of expr in input
                    auto it = find(input.begin(), input.end(), expr);
                    // check if expr is in input
                    if (it >= input.end())
                    {
                        // runtime error
                        throw runtime_error("Error in Function constructor");
                    }
                    // convert it find output to index
                    int index = it - input.begin();
                    AlgEl el({INPUT, index, 0, 0.0});
                    algorithm.push_back(el);
                }
                // check if expr is a Const
                else if (expr->is_const())
                {
                    AlgEl el({CONST, 0, 0, expr->value()});
                    algorithm.push_back(el);
                }
                // check if expr is a Zero
                else if (expr->is_zero())
                {
                    AlgEl el({CONST, 0, 0, 0.0});
                    algorithm.push_back(el);
                }
                // check if expr is a BinaryNode
                else if (expr->is_binary())
                {
                    // find the index of expr1 in ordered_expression
                    auto it = find(ordered_expression.begin(), ordered_expression.end(), expr->dep(0));

                    // check if expr1 is in ordered_expression
                    if (it >= ordered_expression.end())
                    {
                        // runtime error
                        throw runtime_error("Error in Function constructor");
                    }
                    // convert it find output to index
                    int index1 = it - ordered_expression.begin();
                    // find the index of expr2 in ordered_expression
                    it = find(ordered_expression.begin(), ordered_expression.end(), expr->dep(1));
                    // check if expr2 is in ordered_expression
                    if (it >= ordered_expression.end())
                    {
                        // runtime error
                        throw runtime_error("Error in Function constructor");
                    }
                    // convert it find output to index
                    int index2 = it - ordered_expression.begin();
                    if (expr->is_mult())
                    {
                        AlgEl el({MULT, index1, index2, 0.0});
                        algorithm.push_back(el);
                    }
                    else if (expr->is_plus())
                    {
                        AlgEl el({PLUS, index1, index2, 0.0});
                        algorithm.push_back(el);
                    }
                    else
                    {
                        // runtime error
                        throw runtime_error("Error in Function constructor");
                    }
                }
                else
                {
                    // runtime error
                    throw runtime_error("Error in Function constructor");
                }
            }
            // add the output instruction
            for (auto expr : output)
            {
                // find the index of expr in ordered_expression
                auto it = find(ordered_expression.begin(), ordered_expression.end(), expr);
                if (it == ordered_expression.end())
                    throw runtime_error("error in Function constructor");
                int index = it - ordered_expression.begin();
                // find the index of expr in output
                auto it2 = find(output.begin(), output.end(), expr);
                if (it2 == output.end())
                    throw runtime_error("error in Function constructor");
                int index2 = it2 - output.begin();
                AlgEl el({OUTPUT, index, index2, 0.0});
                algorithm.push_back(el);
            }
            output_size = output.size();
        }
        // call the function
        vector<double> Eval(const vector<double> &input)
        {
            vector<double> result(output_size);
            work.resize(0);
            for (auto el : algorithm)
            {
                switch (el.ins)
                {
                case INPUT:
                    work.push_back(input.at(el.arg1));
                    break;
                case OUTPUT:
                    work.push_back(0.0);
                    result.at(el.arg2) = work.at(el.arg1);
                    break;
                case PLUS:
                    work.push_back(work.at(el.arg1) + work.at(el.arg2));
                    break;
                case MULT:
                    work.push_back(work.at(el.arg1) * work.at(el.arg2));
                    break;
                case CONST:
                    work.push_back(el.val);
                    break;
                default:
                    // runtime error
                    throw runtime_error("Error in Function Eval");
                    break;
                }
            }
            return result;
        }
        vector<AlgEl> algorithm;
        vector<double> work;
        vector<Expression> work_e;
        int output_size;
    };

    // matrix stuff
    template <typename T>
    class Matrix
    {
    public:
        Matrix() : Matrix(0, 0){};
        Matrix(const int n_rows, const int n_cols) : n_rows_(n_rows), n_cols_(n_cols)
        {
            data.resize(n_rows * n_cols, T(0.0));
        };
        Matrix(const T &val) : Matrix(1, 1) { data[0] = val; };
        Matrix(const vector<T> &val) : Matrix(val.size(), 1) { data = val; };
        Matrix(const initializer_list<T> &val) : Matrix(val.size(), 1)
        {
            // copy
            copy(val.begin(), val.end(), data.begin());
        };
        T &operator()(int i, int j) { return data[i + n_rows_ * j]; };
        const T &operator()(int i, int j) const { return data[i + n_rows_ * j]; };
        int n_rows() const
        {
            return n_rows_;
        };
        int n_cols() const
        {
            return n_cols_;
        };
        Matrix<T> block(int i, int j, int n_rows, int n_cols) const
        {
            Matrix<T> result(n_rows, n_cols);
            for (int k = 0; k < n_rows; k++)
            {
                for (int l = 0; l < n_cols; l++)
                {
                    result(k, l) = (*this)(i + k, j + l);
                }
            }
            return result;
        };
        Matrix<T> transpose()
        {
            Matrix<T> result(n_cols_, n_rows_);
            for (int i = 0; i < n_rows_; i++)
            {
                for (int j = 0; j < n_cols_; j++)
                {
                    result(j, i) = (*this)(i, j);
                }
            }
            return result;
        }

    protected:
        vector<T> data;
        int n_rows_;
        int n_cols_;
    };

    template <typename T>
    ostream &operator<<(ostream &os, const Matrix<T> &A)
    {
        for (int i = 0; i < A.n_rows(); i++)
        {
            for (int j = 0; j < A.n_cols(); j++)
            {
                os << A(i, j) << " ";
            }
            os << endl;
        }
        return os;
    }
    class Zerom : public Matrix<Expression>
    {
    public:
        Zerom(int n_rows, int n_cols) : Matrix(n_rows, n_cols)
        {
            for (int i = 0; i < n_rows; i++)
            {
                for (int j = 0; j < n_cols; j++)
                {
                    (*this)(i, j) = 0.0;
                }
            }
        };
    };
    template <typename T>
    class Eye : public Matrix<T>
    {
    public:
        Eye(int n_rows) : Matrix<T>(n_rows, n_rows)
        {
            for (int i = 0; i < n_rows; i++)
            {
                for (int j = 0; j < n_rows; j++)
                {
                    if (i == j)
                    {
                        (*this)(i, j) = 1.0;
                    }
                    else
                    {
                        (*this)(i, j) = 0.0;
                    }
                }
            }
        };
    };
    class Constv : public Matrix<Expression>
    {
    public:
        Constv(const double &val) : Matrix(1, 1)
        {
            (*this)(0, 0) = val;
        }
        Constv(const vector<double> &val) : Matrix(val.size(), 1)
        {
            for (int i = 0; i < n_rows(); i++)
            {
                for (int j = 0; j < n_cols(); j++)
                {
                    (*this)(i, j) = val[i + j * n_rows()];
                }
            }
        };
    };

    class SymMatrix : public Matrix<Expression>
    {
    public:
        SymMatrix(const string &name, const int n_rows, const int n_cols) : Matrix<Expression>(n_rows, n_cols)
        {
            for (int i = 0; i < n_rows; i++)
            {
                for (int j = 0; j < n_cols; j++)
                {
                    (*this)(i, j) = Sym(name + string("(") + to_string(i) + string(",") + to_string(j) + string(")"));
                }
            }
        };
    };

    template <typename T>
    Matrix<T> operator*(const Matrix<T> &A, const Matrix<T> &B)
    {
        int n_rows = A.n_rows();
        int n_cols = B.n_cols();
        int n_inner = A.n_cols();
        Zerom C(n_rows, n_cols);
        for (int i = 0; i < n_rows; i++)
        {
            for (int j = 0; j < n_cols; j++)
            {
                for (int k = 0; k < n_inner; k++)
                {
                    C(i, j) = C(i, j) + A(i, k) * B(k, j);
                }
            }
        }
        return C;
    }

    Matrix<Expression> tril(const Matrix<Expression> &A)
    {
        int n_rows = A.n_rows();
        int n_cols = A.n_cols();
        Zerom C(n_rows, n_cols);
        for (int i = 0; i < n_rows; i++)
        {
            for (int j = 0; j < n_cols; j++)
            {
                if (i >= j)
                {
                    C(i, j) = A(i, j);
                }
            }
        }
        return C;
    }

    template <typename T>
    Matrix<T> operator*(const Expression &expr, const Matrix<T> &B)
    {
        int n_rows = B.n_rows();
        int n_cols = B.n_cols();
        Zerom C(n_rows, n_cols);
        for (int i = 0; i < n_rows; i++)
        {
            for (int j = 0; j < n_cols; j++)
            {
                C(i, j) = C(i, j) + expr * B(i, j);
            }
        }
        return C;
    }
    Matrix<Expression> operator-(const Matrix<Expression> &expr1)
    {
        return Const(-1.0) * expr1;
    }

    Matrix<Expression> operator+(const Matrix<Expression> &A, const Matrix<Expression> &B)
    {
        int n_rows = A.n_rows();
        int n_cols = A.n_cols();
        Zerom C(n_rows, n_cols);
        for (int i = 0; i < n_rows; i++)
        {
            for (int j = 0; j < n_cols; j++)
            {
                C(i, j) = A(i, j) + B(i, j);
            }
        }
        return C;
    }
    Matrix<Expression> operator+(const Matrix<Expression> &A, const Constv &B)
    {
        int n_rows = A.n_rows();
        int n_cols = A.n_cols();
        Zerom C(n_rows, n_cols);
        for (int i = 0; i < n_rows; i++)
        {
            for (int j = 0; j < n_cols; j++)
            {
                C(i, j) = A(i, j) + B(i, j);
            }
        }
        return C;
    }
    Matrix<Expression> operator+(const Constv &A, const Matrix<Expression> &B)
    {
        int n_rows = A.n_rows();
        int n_cols = A.n_cols();
        Zerom C(n_rows, n_cols);
        for (int i = 0; i < n_rows; i++)
        {
            for (int j = 0; j < n_cols; j++)
            {
                C(i, j) = A(i, j) + B(i, j);
            }
        }
        return C;
    }

    template <typename T>
    vector<T> vec(const Matrix<T> &mat)
    {
        vector<T> vec;
        for (int j = 0; j < mat.n_cols(); j++)
        {
            for (int i = 0; i < mat.n_rows(); i++)
            {
                vec.push_back(mat(i, j));
            }
        }
        return vec;
    }

    template <typename T>
    vector<T> vec(const vector<Matrix<T>> &mats)
    {
        vector<T> ret;
        for (auto mat : mats)
        {
            auto vecc = vec(mat);
            ret.insert(ret.end(), vecc.begin(), vecc.end());
        }
        return ret;
    }
}
#endif // NODEINCLUDED