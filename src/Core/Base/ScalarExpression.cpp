#include "Core/Base/ScalarExpression.h"
#include <cctype>
#include <cmath>
#include <set>
#include <stdexcept>

namespace Quantum {
namespace {
using Jet=ScalarExpression::Jet;
using Function=std::function<Jet(double,double,double)>;
Jet finite(Jet j){
    if(!std::isfinite(j.value)||!std::isfinite(j.derivative)||!std::isfinite(j.partialSecond))
        throw std::domain_error("analytic expression has a nonfinite value or derivative");
    return j;
}
struct Node {Function evaluate;bool constant=false;};
Node literal(double v){return {[v](double,double,double){return Jet{v,0,0};},true};}
class Parser {
public:
    Parser(const std::string& text,const ScalarExpression::Constants& constants,ScalarExpression::Coordinates coordinates):text_(text),constants_(constants),coordinates_(coordinates){
        if(text_.empty()||text_.size()>4096)fail("formula must contain 1..4096 characters");
        for(const auto& item:constants_){
            if(item.first.empty()||!std::isalpha(static_cast<unsigned char>(item.first[0])))fail("invalid constant name");
            for(char c:item.first)if(!std::isalnum(static_cast<unsigned char>(c))&&c!='_')fail("invalid constant name");
            if(reserved().count(item.first)||variable(item.first)||!std::isfinite(item.second))fail("reserved constant name or nonfinite constant");
        }
    }
    Function parse(){auto n=sum();space();if(at_!=text_.size())fail("unexpected token");return n.evaluate;}
private:
    const std::string& text_;
    const ScalarExpression::Constants& constants_;
    ScalarExpression::Coordinates coordinates_;
    std::size_t at_=0;
    int depth_=0;
    static const std::set<std::string>& reserved(){static const std::set<std::string> names{"s","pi","sin","cos","tan","exp","log","sqrt"};return names;}
    int variable(const std::string& name) const {
        using C=ScalarExpression::Coordinates;
        if(coordinates_==C::Parameter)return name=="s"?1:0;
        if(coordinates_==C::Time)return name=="t"?3:0;
        if((coordinates_==C::PlanarTime || coordinates_==C::MeridionalTime) && name=="t")return 3;
        if(coordinates_==C::Planar || coordinates_==C::PlanarTime)return name=="x"?1:name=="y"?2:0;
        return name=="r"?1:name=="z"?2:0;
    }
    [[noreturn]] void fail(const std::string& message) const {
        throw std::invalid_argument("analytic expression at column "+std::to_string(at_+1)+": "+message);
    }
    void space(){while(at_<text_.size()&&std::isspace(static_cast<unsigned char>(text_[at_])))++at_;}
    bool take(char c){space();if(at_<text_.size()&&text_[at_]==c){++at_;return true;}return false;}
    Node compose(Function f,bool constant){
        if(constant)return literal(finite(f(0,0,0)).value);
        return {std::move(f),false};
    }
    Node sum(){
        auto a=product();
        for(;;){char op=take('+')?'+':take('-')?'-':0;if(!op)return a;auto b=product();
            a=compose([a,b,op](double s,double t,double time){auto x=a.evaluate(s,t,time),y=b.evaluate(s,t,time);const double sign=op=='+'?1.:-1.;
                return finite({x.value+sign*y.value,x.derivative+sign*y.derivative,x.partialSecond+sign*y.partialSecond});},a.constant&&b.constant);
        }
    }
    Node product(){
        auto a=unary();
        for(;;){char op=take('*')?'*':take('/')?'/':0;if(!op)return a;auto b=unary();
            a=compose([a,b,op](double s,double t,double time){auto x=a.evaluate(s,t,time),y=b.evaluate(s,t,time);
                if(op=='*')return finite({x.value*y.value,x.derivative*y.value+x.value*y.derivative,x.partialSecond*y.value+x.value*y.partialSecond});
                if(y.value==0)throw std::domain_error("division by zero in analytic expression");
                const double value=x.value/y.value;
                return finite({value,(x.derivative-value*y.derivative)/y.value,(x.partialSecond-value*y.partialSecond)/y.value});},a.constant&&b.constant);
        }
    }
    Node unary(){
        if(++depth_>64)fail("formula nesting exceeds 64");
        Node result;
        if(take('+'))result=unary();
        else if(take('-')){auto a=unary();result=compose([a](double s,double t,double time){auto x=a.evaluate(s,t,time);return Jet{-x.value,-x.derivative,-x.partialSecond};},a.constant);}
        else result=power();
        --depth_;return result;
    }
    Node power(){
        auto a=primary();if(!take('^'))return a;
        auto b=unary();if(!b.constant)fail("power exponent must be constant");
        const double exponent=b.evaluate(0,0,0).value;
        return compose([a,exponent](double s,double t,double time){auto x=a.evaluate(s,t,time);
            if(exponent==0)return Jet{1,0};
            if(exponent==1)return x;
            const double slope=exponent*std::pow(x.value,exponent-1);
            return finite({std::pow(x.value,exponent),slope*x.derivative,slope*x.partialSecond});},a.constant);
    }
    Node primary(){
        if(take('(')){auto n=sum();if(!take(')'))fail("expected ')'");return n;}
        space();if(at_==text_.size())fail("expected number, variable or function");
        if(std::isdigit(static_cast<unsigned char>(text_[at_]))||text_[at_]=='.'){
            std::size_t used=0;double value=0;
            try{value=std::stod(text_.substr(at_),&used);}catch(const std::exception&){fail("invalid number");}
            if(!used||!std::isfinite(value))fail("number must be finite");at_+=used;return literal(value);
        }
        const auto begin=at_;
        while(at_<text_.size()&&(std::isalnum(static_cast<unsigned char>(text_[at_]))||text_[at_]=='_'))++at_;
        if(at_==begin)fail("expected identifier");
        const auto name=text_.substr(begin,at_-begin);
        const int coordinate=variable(name);
        if(coordinate)return {[coordinate](double s,double t,double time){return coordinate==1?Jet{s,1,0}:coordinate==2?Jet{t,0,1}:Jet{time,0,0};},false};
        if(name=="pi")return literal(std::acos(-1.));
        const auto found=constants_.find(name);if(found!=constants_.end())return literal(found->second);
        if(name=="s"||!reserved().count(name))fail("unknown identifier '"+name+"'");
        if(!take('('))fail("function requires parentheses");auto a=sum();if(!take(')'))fail("expected ')'");
        return compose([a,name](double s,double t,double time){auto x=a.evaluate(s,t,time);double value=0,slope=0;
            if(name=="sin"){value=std::sin(x.value);slope=std::cos(x.value);}
            else if(name=="cos"){value=std::cos(x.value);slope=-std::sin(x.value);}
            else if(name=="tan"){value=std::tan(x.value);slope=1+value*value;}
            else if(name=="exp"){value=std::exp(x.value);slope=value;}
            else if(name=="log"){
                if(!(x.value>0))throw std::domain_error("log requires a positive argument in analytic expression");
                value=std::log(x.value);slope=1/x.value;
            }else {
                if(!(x.value>0))throw std::domain_error("sqrt requires a positive argument for an analytic derivative");
                value=std::sqrt(x.value);slope=.5/value;
            }
            return finite({value,slope*x.derivative,slope*x.partialSecond});},a.constant);
    }
};
}
ScalarExpression::ScalarExpression(const std::string& formula,const Constants& constants,Coordinates coordinates):
    coordinates_(coordinates),evaluate_(Parser(formula,constants,coordinates).parse()){}
ScalarExpression::Jet ScalarExpression::sample(double s) const {
    if(coordinates_!=Coordinates::Parameter && coordinates_!=Coordinates::Time)throw std::invalid_argument("spatial analytic expression needs two coordinates");
    if(coordinates_==Coordinates::Time)return sample(0,0,s);
    return sample(s,0);
}
ScalarExpression::Jet ScalarExpression::sample(double first,double second) const {
    return sample(first,second,0);
}
ScalarExpression::Jet ScalarExpression::sample(double first,double second,double time) const {
    if(!std::isfinite(first)||!std::isfinite(second)||!std::isfinite(time))throw std::invalid_argument("analytic coordinates must be finite");
    return finite(evaluate_(first,second,time));
}
}
