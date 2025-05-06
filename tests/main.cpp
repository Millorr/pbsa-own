#include <iostream>

#include "../Prog1_1Simulation.hpp"
#include "../Prog1_2Simulation.hpp"
#include "../Prog1_3Simulation.hpp"

namespace test_1_1
{
	double potential_energy(Eigen::Vector3d r, double g)
	{
		double const m = 1.0;
		return m * g * r.norm();
	}

	double kinetic_energy(Eigen::Vector3d r, Eigen::Vector3d r_prev, double dt)
	{
		double const m = 1.0;
		auto v2 = (r - r_prev).squaredNorm() / (dt * dt);
		return 0.5 * m * v2;
	}

	int test_explicit(Prog1_1Simulation & simulation)
	{
		std::puts("Running test_1_1::test_explicit...");
		Eigen::Vector3d const r0 = Eigen::Vector3d::UnitX();
		Eigen::Vector3d const dr0 = Eigen::Vector3d::UnitY();
		double const g = 1.0;
		double const dt = 1.0 / 60.0;

		auto p0 = potential_energy(r0, g);
		auto k0 = kinetic_energy(r0, r0 - dt * dr0, dt);

		auto reset = [&] {
			simulation.reset(
				Prog1_1Simulation::Integration::Explicit,
				r0, dr0, g, dt
			);
		};

		{
			reset();
			auto r1 = simulation.step();
			auto r1_ref = Eigen::Vector3d{1.0, 0.01666666666666666666666666666667, 0.0};

			// test single step
			if(!r1.isApprox(r1_ref, 1e-6))
			{
				std::cerr << "ERROR: first step does not match the expected result. " << r1.transpose() << " != " << r1_ref << " (epsilon 1e-6)\n";
				return 1;
			}

			// test that reset is deterministic
			reset();
			auto r1b = simulation.step();

			if(r1 != r1b)
			{
				std::cerr << "ERROR: reset is not deterministic. " << r1.transpose() << " != " << r1b.transpose() << '\n';
				return 1;
			}
		}

		// test one full rotation
		{
			reset();
			auto r_prev = r0;
			auto nq = 0;
			for(int i = 1; i <= 500 && nq < 4; ++i)
			{
				auto r = simulation.step();

				if(std::fabs(r.z()) > 1e-10)
				{
					std::fputs("ERROR: object left the expected plane of motion\n", stderr);
					return 1;
				}

				auto p = potential_energy(r, g);
				auto k = kinetic_energy(r, r_prev, dt);

				if(r.x() <= 0 && r_prev.x() >= 0)
				{
					if(std::abs(i - 95) > 1)
					{
						std::fprintf(stderr, "ERROR: expected object to pass quadrant 0 in 95 iterations (+/- 1), not %d iterations\n", i);
						return 1;
					}
					if(std::fabs(p - 1.02693582412702189) > 1e-3 || std::fabs(k - 0.49957798715595547) > 1e-3)
					{
						std::fprintf(stderr, "ERROR: mismatched potential or kinetic energy after quadrant 0. e_p = %1.17f, e_k = %1.17f, expected e_p = 1.02693582412702189, e_k = 0.49957798715595547 (epsilon 1e-3)\n", p, k);
						return 1;
					}
					++nq;
				}
				else if(r.y() <= 0 && r_prev.y() >= 0)
				{
					if(std::abs(i - 198) > 2)
					{
						std::fprintf(stderr, "ERROR: expected object to pass quadrant 1 in 198 iterations (+/- 2), not %d iterations\n", i);
						return 1;
					}
					if(std::fabs(p - 1.10637435797107564) > 1e-3 || std::fabs(k - 0.45255832683626201) > 1e-3)
					{
						std::fprintf(stderr, "ERROR: mismatched potential or kinetic energy after quadrant 1. e_p = %1.17f, e_k = %1.17f, expected e_p = 1.10637435797107564, e_k = 0.45255832683626201 (epsilon 1e-3)\n", p, k);
						return 1;
					}
					++nq;
				}
				else if(r.x() >= 0 && r_prev.x() <= 0)
				{
					if(std::abs(i - 315) > 3)
					{
						std::fprintf(stderr, "ERROR: expected object to pass quadrant 2 in 315 iterations (+/- 3), not %d iterations\n", i);
						return 1;
					}
					if(std::fabs(p - 1.18584592623053009) > 1e-3 || std::fabs(k - 0.40998074089301895) > 1e-3)
					{
						std::fprintf(stderr, "ERROR: mismatched potential or kinetic energy after quadrant 2. e_p = %1.17f, e_k = %1.17f, expected e_p = 1.18584592623053009, e_k = 0.40998074089301895 (epsilon 1e-3)\n", p, k);
						return 1;
					}
					++nq;
				}
				else if(r.y() >= 0 && r_prev.y() <= 0 && i != 1)
				{
					if(std::abs(i - 440) > 4)
					{
						std::fprintf(stderr, "ERROR: expected object to pass quadrant 3 in 440 iterations (+/- 4), not %d iterations\n", i);
						return 1;
					}
					if(std::fabs(p - 1.19884531983095699) > 1e-3 || std::fabs(k - 0.41729636317557300) > 1e-3)
					{
						std::fprintf(stderr, "ERROR: mismatched potential or kinetic energy after quadrant 3. e_p = %1.17f, e_k = %1.17f, expected e_p = 1.19884531983095699, e_k = 0.41729636317557300 (epsilon 1e-3)\n", p, k);
						return 1;
					}
					++nq;
				}

				r_prev = r;
			}

			if(nq != 4)
			{
				std::fprintf(stderr, "ERROR: object did not complete a full rotation within 500 iterations\n");
				return 1;
			}
		}

		std::puts("test_1_1::test_explicit passed");

		return 0;
	}

	int test_implicit(Prog1_1Simulation & simulation)
	{
		std::puts("Running test_1_1::test_implicit...");
		Eigen::Vector3d const r0 = Eigen::Vector3d::UnitX();
		Eigen::Vector3d const dr0 = Eigen::Vector3d::UnitY();
		double const g = 1.0;
		double const dt = 1.0 / 60.0;

		auto p0 = potential_energy(r0, g);
		auto k0 = kinetic_energy(r0, r0 - dt * dr0, dt);

		auto reset = [&] {
			simulation.reset(
				Prog1_1Simulation::Integration::Implicit,
				r0, dr0, g, dt
			);
		};

		{
			reset();
			auto r1 = simulation.step();
			auto r1_ref = Eigen::Vector3d{0.99972206781545303, 0.01666203832268814, 0.0};

			// test single step
			if(!r1.isApprox(r1_ref, 1e-6))
			{
				std::cerr << "ERROR: first step does not match the expected result. " << r1.transpose() << " != " << r1_ref << " (epsilon 1e-6)\n";
				return 1;
			}

			// test that reset is deterministic
			reset();
			auto r1b = simulation.step();

			if(r1 != r1b)
			{
				std::cerr << "ERROR: reset is not deterministic. " << r1.transpose() << " != " << r1b.transpose() << '\n';
				return 1;
			}
		}

		// test one full rotation
		{
			reset();
			auto r_prev = r0;
			auto nq = 0;
			for(int i = 1; i < 500 && nq < 4; ++i)
			{
				auto r = simulation.step();

				if(std::fabs(r.z()) > 1e-10)
				{
					std::fputs("ERROR: object left the expected plane of motion\n", stderr);
					return 1;
				}

				auto p = potential_energy(r, g);
				auto k = kinetic_energy(r, r_prev, dt);

				if(r.x() <= 0 && r_prev.x() >= 0)
				{
					if(std::abs(i - 94) > 1)
					{
						std::fprintf(stderr, "ERROR: expected object to pass quadrant 0 in 94 iterations (+/- 1), not %d iterations\n", i);
						return 1;
					}
					if(std::fabs(p - 0.97144383724162975) > 1e-3 || std::fabs(k - 0.50237628785036170) > 1e-3)
					{
						std::fprintf(stderr, "ERROR: mismatched potential or kinetic energy after quadrant 0. e_p = %1.17f, e_k = %1.17f, expected e_p = 0.97144383724162975, e_k = 0.50237628785036170 (epsilon 1e-3)\n", p, k);
						return 1;
					}
					++nq;
				}
				else if(r.y() <= 0 && r_prev.y() >= 0)
				{
					if(std::abs(i - 180) > 2)
					{
						std::fprintf(stderr, "ERROR: expected object to pass quadrant 1 in 180 iterations (+/- 2), not %d iterations\n", i);
						return 1;
					}
					if(std::fabs(p - 0.89504346511766952) > 1e-3 || std::fabs(k - 0.55948249798938521) > 1e-3)
					{
						std::fprintf(stderr, "ERROR: mismatched potential or kinetic energy after quadrant 1. e_p = %1.17f, e_k = %1.17f, expected e_p = 0.89504346511766952, e_k = 0.55948249798938521 (epsilon 1e-3)\n", p, k);
						return 1;
					}
					++nq;
				}
				else if(r.x() >= 0 && r_prev.x() <= 0)
				{
					if(std::abs(i - 254) > 3)
					{
						std::fprintf(stderr, "ERROR: expected object to pass quadrant 2 in 254 iterations (+/- 3), not %d iterations\n", i);
						return 1;
					}
					if(std::fabs(p - 0.82003476485104987) > 1e-3 || std::fabs(k - 0.62307963510258690) > 1e-3)
					{
						std::fprintf(stderr, "ERROR: mismatched potential or kinetic energy after quadrant 2. e_p = %1.17f, e_k = %1.17f, expected e_p = 0.82003476485104987, e_k = 0.62307963510258690 (epsilon 1e-3)\n", p, k);
						return 1;
					}
					++nq;
				}
				else if(r.y() >= 0 && r_prev.y() <= 0 && i != 1)
				{
					if(std::abs(i - 321) > 4)
					{
						std::fprintf(stderr, "ERROR: expected object to pass quadrant 3 in 321 iterations (+/- 4), not %d iterations\n", i);
						return 1;
					}
					if(std::fabs(p - 0.77873216435401249) > 1e-3 || std::fabs(k - 0.64153137074583644) > 1e-3)
					{
						std::fprintf(stderr, "ERROR: mismatched potential or kinetic energy after quadrant 3. e_p = %1.17f, e_k = %1.17f, expected e_p = 0.77873216435401249, e_k = 0.64153137074583644 (epsilon 1e-3)\n", p, k);
						return 1;
					}
					++nq;
				}

				r_prev = r;
			}

			if(nq != 4)
			{
				std::fprintf(stderr, "ERROR: object did not complete a full rotation within 500 iterations\n");
				return 1;
			}
		}

		std::puts("test_1_1::test_implicit passed");

		return 0;
	}

	int test_verlet(Prog1_1Simulation & simulation)
	{
		std::puts("Running test_1_1::test_verlet...");
		Eigen::Vector3d const r0 = Eigen::Vector3d::UnitX();
		Eigen::Vector3d const dr0 = Eigen::Vector3d::UnitY();
		double const g = 1.0;
		double const dt = 1.0 / 60.0;

		auto p0 = potential_energy(r0, g);
		auto k0 = kinetic_energy(r0, r0 - dt * dr0, dt);

		auto reset = [&] {
			simulation.reset(
				Prog1_1Simulation::Integration::Verlet,
				r0, dr0, g, dt
			);
		};

		{
			reset();
			auto r1 = simulation.step();
			auto r1_ref = Eigen::Vector3d{0.99986111111111109, 0.01666666666666667, 0.0};

			// test single step
			if(!r1.isApprox(r1_ref, 1e-6))
			{
				std::cerr << "ERROR: first step does not match the expected result. " << r1.transpose() << " != " << r1_ref << " (epsilon 1e-6)\n";
				return 1;
			}

			// test that reset is deterministic
			reset();
			auto r1b = simulation.step();

			if(r1 != r1b)
			{
				std::cerr << "ERROR: reset is not deterministic. " << r1.transpose() << " != " << r1b.transpose() << '\n';
				return 1;
			}
		}

		// test one full rotation
		{
			reset();
			auto r_prev = r0;
			auto nq = 0;
			for(int i = 1; i < 500 && nq < 4; ++i)
			{
				auto r = simulation.step();

				if(std::fabs(r.z()) > 1e-10)
				{
					std::fputs("ERROR: object left the expected plane of motion\n", stderr);
					return 1;
				}

				auto p = potential_energy(r, g);
				auto k = kinetic_energy(r, r_prev, dt);

				if(r.x() <= 0 && r_prev.x() >= 0)
				{
					if(std::abs(i - 95) > 1)
					{
						std::fprintf(stderr, "ERROR: expected object to pass quadrant 0 in 95 iterations (+/- 1), not %d iterations\n", i);
						return 1;
					}
					if(std::fabs(p - 1) > 1e-3 || std::fabs(k - 0.5) > 1e-3)
					{
						std::fprintf(stderr, "ERROR: mismatched potential or kinetic energy after quadrant 0. e_p = %1.17f, e_k = %1.17f, expected e_p = 1, e_k = 0.5 (epsilon 1e-3)\n", p, k);
						return 1;
					}
					++nq;
				}
				else if(r.y() <= 0 && r_prev.y() >= 0)
				{
					if(std::abs(i - 189) > 2)
					{
						std::fprintf(stderr, "ERROR: expected object to pass quadrant 1 in 189 iterations (+/- 2), not %d iterations\n", i);
						return 1;
					}
					if(std::fabs(p - 1) > 1e-3 || std::fabs(k - 0.5) > 1e-3)
					{
						std::fprintf(stderr, "ERROR: mismatched potential or kinetic energy after quadrant 1. e_p = %1.17f, e_k = %1.17f, expected e_p = 1, e_k = 0.5 (epsilon 1e-3)\n", p, k);
						return 1;
					}
					++nq;
				}
				else if(r.x() >= 0 && r_prev.x() <= 0)
				{
					if(std::abs(i - 283) > 3)
					{
						std::fprintf(stderr, "ERROR: expected object to pass quadrant 2 in 283 iterations (+/- 3), not %d iterations\n", i);
						return 1;
					}
					if(std::fabs(p - 1) > 1e-3 || std::fabs(k - 0.5) > 1e-3)
					{
						std::fprintf(stderr, "ERROR: mismatched potential or kinetic energy after quadrant 2. e_p = %1.17f, e_k = %1.17f, expected e_p = 1, e_k = 0.5 (epsilon 1e-3)\n", p, k);
						return 1;
					}
					++nq;
				}
				else if(r.y() >= 0 && r_prev.y() <= 0 && i != 1)
				{
					if(std::abs(i - 378) > 4)
					{
						std::fprintf(stderr, "ERROR: expected object to pass quadrant 3 in 378 iterations (+/- 4), not %d iterations\n", i);
						return 1;
					}
					if(std::fabs(p - 1) > 1e-3 || std::fabs(k - 0.5) > 1e-3)
					{
						std::fprintf(stderr, "ERROR: mismatched potential or kinetic energy after quadrant 3. e_p = %1.17f, e_k = %1.17f, expected e_p = 1, e_k = 0.5 (epsilon 1e-3)\n", p, k);
						return 1;
					}
					++nq;
				}

				r_prev = r;
			}

			if(nq != 4)
			{
				std::fprintf(stderr, "ERROR: object did not complete a full rotation within 500 iterations\n");
				return 1;
			}
		}

		std::puts("test_1_1::test_verlet passed");

		return 0;
	}
}

namespace test_1_2
{
	int test_two_body(Prog1_2Simulation & simulation)
	{
		std::puts("Running test_1_2::test_two_body...");
		std::vector<std::tuple<Eigen::Vector3d, Eigen::Vector3d, double>> const bodies{
			std::make_tuple(Eigen::Vector3d{ 0.5, 0, 0}, Eigen::Vector3d{0,  0.7071067812, 0}, 1.0),
			std::make_tuple(Eigen::Vector3d{-0.5, 0, 0}, Eigen::Vector3d{0, -0.7071067812, 0}, 1.0)
		};
		auto const G = 1;
		auto const dt = 0.05;

		auto reset = [&] {
			simulation.reset(bodies, G, dt);
		};

		{
			reset();
			// test single step
			auto xs1 = simulation.step();
			if(xs1.size() != 2)
			{
				std::fprintf(stderr, "ERROR: expected 2 positions, got %zu\n", xs1.size());
				return 1;
			}
			auto x0_ref = Eigen::Vector3d{0.49875, 0.03535533906, 0.0};
			auto x1_ref = Eigen::Vector3d{-0.49875, -0.03535533906, 0.0};
			Eigen::Matrix<double, 3, 2> xs, xs_ref;
			xs << xs1[0], xs1[1];
			xs_ref << x0_ref, x1_ref;
			if(!xs.isApprox(xs_ref, 1e-3))
			{
				std::cerr
					<< "ERROR: first step does not match the expected result.\n"
					<< xs
					<< "\n!=\n"
					<< xs_ref
					<< "\n(epsilon 1e-3)\n";
				return 1;
			}

			// test that reset is deterministic
			reset();
			auto xs1b = simulation.step();
			if(xs1 != xs1b)
			{
				Eigen::Matrix<double, 3, 2> xsb;
				xsb << xs1b[0], xs1b[1];
				std::cerr
					<< "ERROR: reset is not deterministic.\n"
					<< xs
					<< "\n!=\n"
					<< xsb
					<< '\n';
				return 1;
			}
		}

		// test a full revolution
		{
			reset();
			auto x0_prev = std::get<0>(bodies[0]);

			bool completed = false;
			for(int i = 1; i <= 500 && !completed; ++i)
			{
				auto const & xs = simulation.step();
				if(xs.size() != 2)
				{
					std::fprintf(stderr, "ERROR: expected 2 positions, got %zu\n", xs.size());
					return 1;
				}

				auto x0 = xs[0];
				auto x1 = xs[1];

				if(std::fabs(x0.z()) > 1e-5 || std::fabs(x1.z()) > 1e-5)
				{
					std::fputs("ERROR: objects left the expected plane of motion\n", stderr);
					return 1;
				}

				if(!x0.isApprox(-x1, 1e-5))
				{
					std::fputs("ERROR: objects have left their expected center of mass\n", stderr);
					return 1;
				}

				if(x0.y() >= 0 && x0_prev.y() <= 0 && i != 1)
				{
					if(std::abs(i - 90) > 2)
					{
						std::cout << x0_prev.transpose() << ", " << x0.transpose() << '\n';
						std::fprintf(stderr, "ERROR: expected objects to complete a mutual orbit in 90 iterations (+/- 2), not %d iterations\n", i);
						return 1;
					}
					completed = true;
				}

				x0_prev = x0;
			}

			if(!completed)
			{
				std::fprintf(stderr, "ERROR: objects did not complete a full mutual orbit within 500 iterations\n");
				return 1;
			}
		}

		std::puts("test_1_2::test_two_body passed");

		return 0;
	}

	int test_three_body(Prog1_2Simulation & simulation)
	{
		std::puts("Running test_1_2::test_three_body...");
		std::vector<std::tuple<Eigen::Vector3d, Eigen::Vector3d, double>> const bodies{
			std::make_tuple(Eigen::Vector3d{ 0.5, 0, 0}, Eigen::Vector3d{0,  0.7071067812, -5e-5}, 1e4),
			std::make_tuple(Eigen::Vector3d{-0.5, 0, 0}, Eigen::Vector3d{0, -0.7071067812, -5e-5}, 1e4),
			std::make_tuple(Eigen::Vector3d{ 2.5, 0, 0}, Eigen::Vector3d{0, 0, 1}, 1.0),
		};
		auto G = 1e-4;
		auto const dt = 0.2;

		auto reset = [&] {
			simulation.reset(bodies, G, dt);
		};

		{
			reset();
			// test single step
			auto xs1 = simulation.step();
			if(xs1.size() != 3)
			{
				std::fprintf(stderr, "ERROR: expected 3 positions, got %zu\n", xs1.size());
				return 1;
			}
			auto x0_ref = Eigen::Vector3d{0.4800004999999999, 0.14142135624, -1e-5};
			auto x1_ref = Eigen::Vector3d{-0.4799997777777777, -0.14142135624, -1e-5};
			auto x2_ref = Eigen::Vector3d{2.4927777777777775, 0.0, 0.2};
			Eigen::Matrix3d xs, xs_ref;
			xs << xs1[0], xs1[1], xs1[2];
			xs_ref << x0_ref, x1_ref, x2_ref;
			if(!xs.isApprox(xs_ref, 1e-3))
			{
				std::cerr
					<< "ERROR: first step does not match the expected result.\n"
					<< xs
					<< "\n!=\n"
					<< xs_ref
					<< "\n(epsilon 1e-3)\n";
				return 1;
			}

			// test that reset is deterministic
			reset();
			auto xs1b = simulation.step();

			if(xs1 != xs1b)
			{
				Eigen::Matrix3d xsb;
				xsb << xs1b[0], xs1b[1], xs1b[2];
				std::cerr
					<< "ERROR: reset is not deterministic.\n"
					<< xs
					<< "\n!=\n"
					<< xsb
					<< '\n';
				return 1;
			}
		}

		// test a full revolution
		{
			reset();
			auto x0_prev = std::get<0>(bodies[0]);
			auto x2_prev = std::get<0>(bodies[2]);

			bool completed_inner = false;
			bool completed_outer = false;
			for(int i = 1; i <= 500 && !(completed_inner && completed_outer); ++i)
			{
				auto const & xs = simulation.step();
				if(xs.size() != 3)
				{
					std::fprintf(stderr, "ERROR: expected 3 positions, got %zu\n", xs.size());
					return 1;
				}

				auto x0 = xs[0];
				auto x1 = xs[1];
				auto x2 = xs[2];

				// allow fairly large errors, as the objects do not follow a strictly planar orbit
				if(std::fabs(x0.z()) > 1e-2 || std::fabs(x1.z()) > 1e-2 || std::fabs(x2.y()) > 1e-1)
				{
					std::fputs("ERROR: objects left the expected plane of motion\n", stderr);
					std::cout << x0.transpose() << '\n';
					std::cout << x1.transpose() << '\n';
					std::cout << x2.transpose() << '\n';
					return 1;
				}

				if(!x0.isApprox(-x1, 1e-2))
				{
					std::fputs("ERROR: objects have left their expected center of mass\n", stderr);
					std::cout << x0.transpose() << '\n';
					std::cout << x1.transpose() << '\n';
					return 1;
				}

				if(x0.y() >= 0 && x0_prev.y() <= 0 && i != 1 && !completed_inner)
				{
					if(std::abs(i - 23) > 3)
					{
						std::cout << x0_prev.transpose() << ", " << x0.transpose() << '\n';
						std::fprintf(stderr, "ERROR: expected inner objects to complete a mutual orbit in 23 iterations (+/- 3), not %d iterations\n", i);
						return 1;
					}
					completed_inner = true;
				}

				if(x2.z() >= 0 && x2_prev.z() <= 0 && i != 1 && !completed_outer)
				{
					if(std::abs(i - 132) > 4)
					{
						std::cout << x2_prev.transpose() << ", " << x2.transpose() << '\n';
						std::fprintf(stderr, "ERROR: expected outer object to complete an orbit in 132 iterations (+/- 4), not %d iterations\n", i);
						return 1;
					}
					completed_outer = true;
				}

				x0_prev = x0;
				x2_prev = x2;
			}

			if(!completed_inner || !completed_outer)
			{
				std::fprintf(stderr, "ERROR: objects did not complete a full mutual orbit within 500 iterations\n");
				return 1;
			}
		}

		std::puts("test_1_2::test_three_body passed");
		return 0;
	}
}

namespace test_1_3
{
	constexpr std::size_t GRID_LENGTH = 101;

	int test_heat_sim(Prog1_3Simulation & simulation)
	{
		std::puts("Running test_1_3::test_heat_sim...");

		auto const a = 100;
		auto const dt = 0.01;

		auto reset = [&] {
			simulation.reset(dt, a);
		};

		int i_center = GRID_LENGTH / 2, j_center = GRID_LENGTH / 2;

		// one step
		{
			reset();
			auto u1 = simulation.step();

			// check size of return value
			if(u1.size() != GRID_LENGTH * GRID_LENGTH)
			{
				std::fprintf(stderr, "ERROR: Output size and GRID_LENGTH mismatch. Expected: %zi, not %zi", GRID_LENGTH * GRID_LENGTH, u1.size());
				return 1;
			}

			if(std::fabs(u1[i_center * GRID_LENGTH + j_center] - 0.254050) > 1e-3)
			{
				std::fprintf(stderr, "ERROR: Incorrect heat value at center after first step. Expected: %f, not %f. \n", 0.254050, u1[i_center * GRID_LENGTH + j_center]);
				return 1;
			}

			// test that reset is deterministic
			reset();
			auto u1b = simulation.step();

			if(u1 != u1b)
			{
				std::cerr << "ERROR: reset is not deterministic.\n";
				return 1;
			}
		}



		// simulate 500 steps
		{
			reset();
			auto u = simulation.step();
			double sumPrev = 0.;
			for(int i = 0; i < 500; ++i) 
			{
				sumPrev = u.sum();
				u = simulation.step();
				// check total heat
				if(!(u.sum() <= sumPrev))
				{
					if(std::fabs(u.sum() - sumPrev) < 1e-6)
						continue;
					std::fprintf(stderr, "ERROR: Energy gain after %i steps. Total heat now: %f was %f \n", i, u.sum(), sumPrev);
					return 1;
				}
			}

			// check dirichlet boundary
			double sumBoundary = 0.;

			for(std::size_t i = 0; i < GRID_LENGTH; ++i)
			{
				int j_top = 0;
				int j_bottom = GRID_LENGTH - 1;
				sumBoundary += u[i * GRID_LENGTH + j_top]; // Top
				sumBoundary += u[i * GRID_LENGTH + j_bottom]; // Bottom
			}

			for(std::size_t j = 1; j < GRID_LENGTH - 1; ++j)
			{
				int i_left = 0;
				int i_right = GRID_LENGTH - 1;
				sumBoundary += u[i_left * GRID_LENGTH + j]; // Left
				sumBoundary += u[i_right * GRID_LENGTH + j]; // Right
			}

			if(std::fabs(sumBoundary) > 0)
			{
				std::fprintf(stderr, "ERROR: Dirichlet boundary not zero. Total heat on boundary = %f", sumBoundary);
				return 1;
			}

			// check total heat value
			if(std::fabs(u.sum() - 0.595681) > 1e-3)
			{
				std::fprintf(stderr, "ERROR:. Total heat after 500 simulation steps: %f expected %f \n", u.sum(), 0.595681);
				return 1;
			}

		}

		std::puts("test_1_3::test_heat_sim passed");
		return 0;
	}
}

int main(int, char **)
{
	int err = 0;
	int shift = 0;

	{
		Prog1_1Simulation simulation;
		err |= test_1_1::test_explicit(simulation) << (shift++);
		err |= test_1_1::test_implicit(simulation) << (shift++);
		err |= test_1_1::test_verlet(simulation)   << (shift++);
	}

	{
		Prog1_2Simulation simulation;
		err |= test_1_2::test_two_body(simulation)  << (shift++);
		err |= test_1_2::test_three_body(simulation) << (shift++);
	}

	{
		Prog1_3Simulation simulation;
		err |= test_1_3::test_heat_sim(simulation) << (shift++);
	}

	if(!err)
		std::puts("All tests passed");

	return err;
}
