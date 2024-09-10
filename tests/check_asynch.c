#include <stdlib.h>
#include <check.h>

#include <date_manip.h>

START_TEST (test_date_manip_days_in_month)
{
    ck_assert( days_in_month(0,1995) == 31 );
    ck_assert( days_in_month(1,1995) == 28);
    ck_assert( days_in_month(2,1995) == 31 );
    ck_assert( days_in_month(3,1995) == 30);

    ck_assert( days_in_month(1,2000) == 29);
}
END_TEST


Suite * asynch_suite(void)
{
    Suite *s;
    TCase *tc_date_manip;

    s = suite_create("Asynch");

    /* Date manip test case */
    tc_date_manip = tcase_create("Date manip ");

    tcase_add_test(tc_date_manip , test_date_manip_days_in_month);
    suite_add_tcase(s, tc_date_manip );

    return s;
}


int main(void)
{
    int number_failed;
    Suite *s;
    SRunner *sr;

    s = asynch_suite();
    sr = srunner_create(s);

    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}