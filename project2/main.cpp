/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*!
Copyright (C) 2005, 2006, 2007, 2009 StatPro Italia srl

This file is part of QuantLib, a free-software/open-source library
for financial quantitative analysts and developers - http://quantlib.org/

QuantLib is free software: you can redistribute it and/or modify it
under the terms of the QuantLib license.  You should have received a
copy of the license along with this program; if not, please email
<quantlib-dev@lists.sf.net>. The license is also available online at
<http://quantlib.org/license.shtml>.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#include <ql/qldefines.hpp>
#ifdef BOOST_MSVC
#  include <ql/auto_link.hpp>
#endif
#include <ql/instruments/vanillaoption.hpp>
#include <ql/pricingengines/vanilla/binomialengine.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/utilities/dataformatters.hpp>
#include <ql/termstructures/yield/zerocurve.hpp>
#include <ql/termstructures/volatility/equityfx/blackvariancecurve.hpp>
#include <boost/timer.hpp>
#include <iostream>
#include <stdlib.h>
#include <iomanip>

using namespace QuantLib;

#if defined(QL_ENABLE_SESSIONS)
namespace QuantLib {

	Integer sessionId() { return 0; }

}
#endif

inline void print_time(boost::timer& timer) {
	double seconds = timer.elapsed();
	Integer hours = int(seconds / 3600);
	seconds -= hours * 3600;
	Integer minutes = int(seconds / 60);
	seconds -= minutes * 60;
	Integer mseconds = int(1000 * (seconds - int(seconds)));
	std::cout << " \nRun completed in ";
	if (hours > 0)
		std::cout << hours << " h ";
	if (hours > 0 || minutes > 0)
		std::cout << minutes << " m ";
	if (hours > 0 || minutes > 0 || int(seconds) > 0)
		std::cout << int(seconds) << "s ";
	std::cout << std::fixed << mseconds << "ms\n" << std::endl;
	timer.restart();

}

int main(int, char*[]) {

	try {

		boost::timer timer;
		std::cout << std::endl;

		// set up dates
		Calendar calendar = TARGET();
		Date todaysDate(15, May, 1998);
		Date settlementDate(17, May, 1998);
		Settings::instance().evaluationDate() = todaysDate;

		// our options
		Option::Type type(Option::Put);
		Real underlying = 36;
		Real strike = 40;
		Spread dividendYield = 0.00;
		Rate riskFreeRate = 0.06;
		Volatility volatility = 0.20;
		Date maturity(17, May, 1999);
		DayCounter dayCounter = Actual365Fixed();

		std::cout << "Option type = " << type << std::endl;
		std::cout << "Maturity = " << maturity << std::endl;
		std::cout << "Underlying price = " << underlying << std::endl;
		std::cout << "Strike = " << strike << std::endl;
		std::cout << "Risk-free interest rate = " << io::rate(riskFreeRate)
			<< std::endl;
		std::cout << "Dividend yield = " << io::rate(dividendYield)
			<< std::endl;
		std::cout << "Volatility = " << io::volatility(volatility)
			<< std::endl;
		std::cout << std::endl;
		std::string method;
		std::cout << std::endl;

		// write column headings
		Size widths[] = { 35, 14, 14, 14 };
		std::cout << std::setw(widths[0]) << std::left << "Method"
			<< std::setw(widths[1]) << std::left << "European"
			<< std::setw(widths[2]) << std::left << "Bermudan"
			<< std::setw(widths[3]) << std::left << "American"
			<< std::endl;

		std::vector<Date> exerciseDates;
		for (Integer i = 1; i <= 4; i++)
			exerciseDates.push_back(settlementDate + 3 * i*Months);

		boost::shared_ptr<Exercise> europeanExercise(
			new EuropeanExercise(maturity));

		boost::shared_ptr<Exercise> bermudanExercise(
			new BermudanExercise(exerciseDates));

		boost::shared_ptr<Exercise> americanExercise(
			new AmericanExercise(settlementDate,
			maturity));

		Handle<Quote> underlyingH(
			boost::shared_ptr<Quote>(new SimpleQuote(underlying)));

		// bootstrap the yield/dividend/vol curves
		//flat yield curve

		Handle<YieldTermStructure> flatTermStructure(
			boost::shared_ptr<YieldTermStructure>(
			new FlatForward(settlementDate, riskFreeRate, dayCounter)));

		// non constant yield curve
		std::vector<Date> dates(3);
		dates[0] = settlementDate;
		dates[1] = settlementDate + Period(1, Years);
		dates[2] = settlementDate + Period(2, Years);
		std::vector<Rate> rates(3);
		rates[0] = 0.05;
		rates[1] = 0.06;
		rates[2] = 0.075;
		Handle<YieldTermStructure> TermStructure(
			boost::shared_ptr<YieldTermStructure>(
			new ZeroCurve(dates, rates, dayCounter)));

		//flat divs curve
		Handle<YieldTermStructure> flatDividendTS(
			boost::shared_ptr<YieldTermStructure>(
			new FlatForward(settlementDate, dividendYield, dayCounter)));

		//flat vol curve
		Handle<BlackVolTermStructure> flatVolTS(
			boost::shared_ptr<BlackVolTermStructure>(
			new BlackConstantVol(settlementDate, calendar, volatility,
			dayCounter)));

		//non constant vol curve
		std::vector<Date> voldates(3);
		voldates[0] = settlementDate + Period(4, Months);
		voldates[1] = settlementDate + Period(8, Months);
		voldates[2] = settlementDate + Period(1, Years);
		std::vector<Volatility> volatilities(3);
		volatilities[0] = 0.025;
		volatilities[1] = 0.036;
		volatilities[2] = 0.0475;
		Handle<BlackVolTermStructure> VolTS(
			boost::shared_ptr<BlackVolTermStructure>(
			new BlackVarianceCurve(settlementDate,voldates, volatilities,
			dayCounter)));

		//payoff and process
		boost::shared_ptr<StrikedTypePayoff> payoff(
			new PlainVanillaPayoff(type, strike));

		boost::shared_ptr<BlackScholesMertonProcess> bsmProcess(
			new BlackScholesMertonProcess(underlyingH, flatDividendTS,
			TermStructure, VolTS)); // use TermStructure, VolTS instead for non constant parameters

		// options
		VanillaOption europeanOption(payoff, europeanExercise);
		VanillaOption bermudanOption(payoff, bermudanExercise);
		VanillaOption americanOption(payoff, americanExercise);

		// Binomial method: Jarrow-Rudd
		Size timeSteps = 3000;

		method = "Binomial Jarrow-Rudd";
		europeanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<JarrowRudd>(bsmProcess, timeSteps)));
		bermudanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<JarrowRudd>(bsmProcess, timeSteps)));
		americanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<JarrowRudd>(bsmProcess, timeSteps)));
		std::cout << std::setw(widths[0]) << std::left << method
			<< std::fixed
			<< std::setw(widths[1]) << std::left << europeanOption.NPV()
			<< std::setw(widths[2]) << std::left << bermudanOption.NPV()
			<< std::setw(widths[3]) << std::left << americanOption.NPV()
			<< std::endl;

		print_time(timer);

		method = "Extended Binomial Jarrow-Rudd";
		europeanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedJarrowRudd_2>(bsmProcess, timeSteps)));
		bermudanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedJarrowRudd_2>(bsmProcess, timeSteps)));
		americanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedJarrowRudd_2>(bsmProcess, timeSteps)));
		std::cout << std::setw(widths[0]) << std::left << method
		<< std::fixed
		<< std::setw(widths[1]) << std::left << europeanOption.NPV()
		<< std::setw(widths[2]) << std::left << bermudanOption.NPV()
		<< std::setw(widths[3]) << std::left << americanOption.NPV()
		<< std::endl;

		print_time(timer);

		method = "Cached Binomial Jarrow-Rudd";
		europeanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedJarrowRudd>(bsmProcess, timeSteps)));
		bermudanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedJarrowRudd>(bsmProcess, timeSteps)));
		americanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedJarrowRudd>(bsmProcess, timeSteps)));
		std::cout << std::setw(widths[0]) << std::left << method
			<< std::fixed
			<< std::setw(widths[1]) << std::left << europeanOption.NPV()
			<< std::setw(widths[2]) << std::left << bermudanOption.NPV()
			<< std::setw(widths[3]) << std::left << americanOption.NPV()
			<< std::endl;

		print_time(timer);

		method = "Binomial Cox-Ross-Rubinstein";
		europeanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<CoxRossRubinstein>(bsmProcess,
			timeSteps)));
		bermudanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<CoxRossRubinstein>(bsmProcess,
			timeSteps)));
		americanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<CoxRossRubinstein>(bsmProcess,
			timeSteps)));
		std::cout << std::setw(widths[0]) << std::left << method
			<< std::fixed
			<< std::setw(widths[1]) << std::left << europeanOption.NPV()
			<< std::setw(widths[2]) << std::left << bermudanOption.NPV()
			<< std::setw(widths[3]) << std::left << americanOption.NPV()
			<< std::endl;

		print_time(timer);

		method = "Extended Binomial Cox-Ross-Rubinstein";
		europeanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedCoxRossRubinstein_2>(bsmProcess,
		timeSteps)));
		bermudanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedCoxRossRubinstein_2>(bsmProcess,
		timeSteps)));
		americanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedCoxRossRubinstein_2>(bsmProcess,
		timeSteps)));
		std::cout << std::setw(widths[0]) << std::left << method
		<< std::fixed
		<< std::setw(widths[1]) << std::left << europeanOption.NPV()
		<< std::setw(widths[2]) << std::left << bermudanOption.NPV()
		<< std::setw(widths[3]) << std::left << americanOption.NPV()
		<< std::endl;

		print_time(timer);

		method = "Cached Binomial Cox-Ross-Rubinstein";
		europeanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedCoxRossRubinstein>(bsmProcess,
			timeSteps)));
		bermudanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedCoxRossRubinstein>(bsmProcess,
			timeSteps)));
		americanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedCoxRossRubinstein>(bsmProcess,
			timeSteps)));
		std::cout << std::setw(widths[0]) << std::left << method
			<< std::fixed
			<< std::setw(widths[1]) << std::left << europeanOption.NPV()
			<< std::setw(widths[2]) << std::left << bermudanOption.NPV()
			<< std::setw(widths[3]) << std::left << americanOption.NPV()
			<< std::endl;

		print_time(timer);

		// Binomial method: Additive equiprobabilities

		method = "Additive equiprobabilities";
		europeanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<AdditiveEQPBinomialTree>(bsmProcess,
			timeSteps)));
		bermudanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<AdditiveEQPBinomialTree>(bsmProcess,
			timeSteps)));
		americanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<AdditiveEQPBinomialTree>(bsmProcess,
			timeSteps)));
		std::cout << std::setw(widths[0]) << std::left << method
			<< std::fixed
			<< std::setw(widths[1]) << std::left << europeanOption.NPV()
			<< std::setw(widths[2]) << std::left << bermudanOption.NPV()
			<< std::setw(widths[3]) << std::left << americanOption.NPV()
			<< std::endl;

		print_time(timer);


		method = "Extended Additive equiprobabilities";
		europeanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedAdditiveEQPBinomialTree_2>(bsmProcess,
		timeSteps)));
		bermudanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedAdditiveEQPBinomialTree_2>(bsmProcess,
		timeSteps)));
		americanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedAdditiveEQPBinomialTree_2>(bsmProcess,
		timeSteps)));
		std::cout << std::setw(widths[0]) << std::left << method
		<< std::fixed
		<< std::setw(widths[1]) << std::left << europeanOption.NPV()
		<< std::setw(widths[2]) << std::left << bermudanOption.NPV()
		<< std::setw(widths[3]) << std::left << americanOption.NPV()
		<< std::endl;

		print_time(timer);

		method = "Cached Additive equiprobabilities";
		europeanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedAdditiveEQPBinomialTree>(bsmProcess,
			timeSteps)));
		bermudanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedAdditiveEQPBinomialTree>(bsmProcess,
			timeSteps)));
		americanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedAdditiveEQPBinomialTree>(bsmProcess,
			timeSteps)));
		std::cout << std::setw(widths[0]) << std::left << method
			<< std::fixed
			<< std::setw(widths[1]) << std::left << europeanOption.NPV()
			<< std::setw(widths[2]) << std::left << bermudanOption.NPV()
			<< std::setw(widths[3]) << std::left << americanOption.NPV()
			<< std::endl;

		print_time(timer);

		// Binomial method: Binomial Trigeorgis

		method = "Binomial Trigeorgis";
		europeanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<Trigeorgis>(bsmProcess, timeSteps)));
		bermudanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<Trigeorgis>(bsmProcess, timeSteps)));
		americanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<Trigeorgis>(bsmProcess, timeSteps)));
		std::cout << std::setw(widths[0]) << std::left << method
			<< std::fixed
			<< std::setw(widths[1]) << std::left << europeanOption.NPV()
			<< std::setw(widths[2]) << std::left << bermudanOption.NPV()
			<< std::setw(widths[3]) << std::left << americanOption.NPV()
			<< std::endl;

		print_time(timer);

		method = "Extended Binomial Trigeorgis";
		europeanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedTrigeorgis_2>(bsmProcess, timeSteps)));
		bermudanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedTrigeorgis_2>(bsmProcess, timeSteps)));
		americanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedTrigeorgis_2>(bsmProcess, timeSteps)));
		std::cout << std::setw(widths[0]) << std::left << method
		<< std::fixed
		<< std::setw(widths[1]) << std::left << europeanOption.NPV()
		<< std::setw(widths[2]) << std::left << bermudanOption.NPV()
		<< std::setw(widths[3]) << std::left << americanOption.NPV()
		<< std::endl;

		print_time(timer);

		method = "Cached Binomial Trigeorgis";
		europeanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedTrigeorgis>(bsmProcess, timeSteps)));
		bermudanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedTrigeorgis>(bsmProcess, timeSteps)));
		americanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedTrigeorgis>(bsmProcess, timeSteps)));
		std::cout << std::setw(widths[0]) << std::left << method
			<< std::fixed
			<< std::setw(widths[1]) << std::left << europeanOption.NPV()
			<< std::setw(widths[2]) << std::left << bermudanOption.NPV()
			<< std::setw(widths[3]) << std::left << americanOption.NPV()
			<< std::endl;

		print_time(timer);

		// Binomial method: Binomial Tian

		method = "Binomial Tian";
		europeanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<Tian>(bsmProcess, timeSteps)));
		bermudanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<Tian>(bsmProcess, timeSteps)));
		americanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<Tian>(bsmProcess, timeSteps)));
		std::cout << std::setw(widths[0]) << std::left << method
			<< std::fixed
			<< std::setw(widths[1]) << std::left << europeanOption.NPV()
			<< std::setw(widths[2]) << std::left << bermudanOption.NPV()
			<< std::setw(widths[3]) << std::left << americanOption.NPV()
			<< std::endl;

		print_time(timer);

		method = "Extended Binomial Tian";
		europeanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedTian_2>(bsmProcess, timeSteps)));
		bermudanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedTian_2>(bsmProcess, timeSteps)));
		americanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedTian_2>(bsmProcess, timeSteps)));
		std::cout << std::setw(widths[0]) << std::left << method
		<< std::fixed
		<< std::setw(widths[1]) << std::left << europeanOption.NPV()
		<< std::setw(widths[2]) << std::left << bermudanOption.NPV()
		<< std::setw(widths[3]) << std::left << americanOption.NPV()
		<< std::endl;

		print_time(timer);

		method = "Cached Binomial Tian";
		europeanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedTian>(bsmProcess, timeSteps)));
		bermudanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedTian>(bsmProcess, timeSteps)));
		americanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedTian>(bsmProcess, timeSteps)));
		std::cout << std::setw(widths[0]) << std::left << method
			<< std::fixed
			<< std::setw(widths[1]) << std::left << europeanOption.NPV()
			<< std::setw(widths[2]) << std::left << bermudanOption.NPV()
			<< std::setw(widths[3]) << std::left << americanOption.NPV()
			<< std::endl;

		print_time(timer);

		// Binomial method: Binomial Leisen-Reimer

		method = "Binomial Leisen-Reimer";
		europeanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<LeisenReimer>(bsmProcess, timeSteps)));
		bermudanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<LeisenReimer>(bsmProcess, timeSteps)));
		americanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<LeisenReimer>(bsmProcess, timeSteps)));
		std::cout << std::setw(widths[0]) << std::left << method
			<< std::fixed
			<< std::setw(widths[1]) << std::left << europeanOption.NPV()
			<< std::setw(widths[2]) << std::left << bermudanOption.NPV()
			<< std::setw(widths[3]) << std::left << americanOption.NPV()
			<< std::endl;

		print_time(timer);

		method = "Extended Binomial Leisen-Reimer";
		europeanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedLeisenReimer_2>(bsmProcess, timeSteps)));
		bermudanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedLeisenReimer_2>(bsmProcess, timeSteps)));
		americanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedLeisenReimer_2>(bsmProcess, timeSteps)));
		std::cout << std::setw(widths[0]) << std::left << method
		<< std::fixed
		<< std::setw(widths[1]) << std::left << europeanOption.NPV()
		<< std::setw(widths[2]) << std::left << bermudanOption.NPV()
		<< std::setw(widths[3]) << std::left << americanOption.NPV()
		<< std::endl;

		print_time(timer);

		method = "Cached Binomial Leisen-Reimer";
		europeanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedLeisenReimer>(bsmProcess, timeSteps)));
		bermudanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedLeisenReimer>(bsmProcess, timeSteps)));
		americanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedLeisenReimer>(bsmProcess, timeSteps)));
		std::cout << std::setw(widths[0]) << std::left << method
			<< std::fixed
			<< std::setw(widths[1]) << std::left << europeanOption.NPV()
			<< std::setw(widths[2]) << std::left << bermudanOption.NPV()
			<< std::setw(widths[3]) << std::left << americanOption.NPV()
			<< std::endl;

		print_time(timer);

		// Binomial method: Binomial Joshi

		method = "Binomial Joshi";
		europeanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<Joshi4>(bsmProcess, timeSteps)));
		bermudanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<Joshi4>(bsmProcess, timeSteps)));
		americanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<Joshi4>(bsmProcess, timeSteps)));
		std::cout << std::setw(widths[0]) << std::left << method
			<< std::fixed
			<< std::setw(widths[1]) << std::left << europeanOption.NPV()
			<< std::setw(widths[2]) << std::left << bermudanOption.NPV()
			<< std::setw(widths[3]) << std::left << americanOption.NPV()
			<< std::endl;

		print_time(timer);

		method = "Extended Binomial Joshi";
		europeanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedJoshi4_2>(bsmProcess, timeSteps)));
		bermudanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedJoshi4_2>(bsmProcess, timeSteps)));
		americanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedJoshi4_2>(bsmProcess, timeSteps)));
		std::cout << std::setw(widths[0]) << std::left << method
		<< std::fixed
		<< std::setw(widths[1]) << std::left << europeanOption.NPV()
		<< std::setw(widths[2]) << std::left << bermudanOption.NPV()
		<< std::setw(widths[3]) << std::left << americanOption.NPV()
		<< std::endl;

		print_time(timer);

		method = "Cached Binomial Joshi";
		europeanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedJoshi4>(bsmProcess, timeSteps)));
		bermudanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedJoshi4>(bsmProcess, timeSteps)));
		americanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<ExtendedJoshi4>(bsmProcess, timeSteps)));
		std::cout << std::setw(widths[0]) << std::left << method
			<< std::fixed
			<< std::setw(widths[1]) << std::left << europeanOption.NPV()
			<< std::setw(widths[2]) << std::left << bermudanOption.NPV()
			<< std::setw(widths[3]) << std::left << americanOption.NPV()
			<< std::endl;

		print_time(timer);

		// End test

		system("pause");
		return 0;

	}
	catch (std::exception& e) {
		std::cerr << e.what() << std::endl;
		return 1;
	}
	catch (...) {
		std::cerr << "unknown error" << std::endl;
		return 1;
	}
}
