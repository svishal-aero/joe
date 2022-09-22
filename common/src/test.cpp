#include <iostream>

enum PrcommActionEnum
{
    REPLACE_DATA=1,
    ADD_DATA,
    MIN_DATA,
    MAX_DATA,
    REPLACE_TRANSLATE_DATA,
    REPLACE_ROTATE_DATA,
    ADD_ROTATE_DATA,
    BITWISE_OR_DATA,
    MIN_NO_PERIODIC_DATA,
    MAX_NO_PERIODIC_DATA,
    ADD_NO_PERIODIC_DATA,
    ADD_TRANSLATE_DATA,
    SUBTRACT_DATA,
    SUBTRACT_ROTATE_DATA,
    NO_CHANGE_DATA
};

enum PrcommTagEnum
{
    UPDATE_R1_TAG=1002,
    UPDATE_R2_TAG,
    UPDATE_R2_REVERSE_TAG,
    EXCHANGE_INT_TAG,
    EXCHANGE_INT_REVERSE_TAG,
    UPDATE_I1_TAG,
    UPDATE_I1_REVERSE_TAG,
    UPDATE_I2_TAG,
    UPDATE_I2_REVERSE_TAG,
    UPDATE_R3_TAG,
    UPDATE_SYMMETRIC_R3_TAG,
    EXCHANGE_PRCOMM_INTV_TAG1,
    EXCHANGE_PRCOMM_INTV_TAG2,
    EXCHANGE_PRCOMM_DOUBLEV_TAG1,
    EXCHANGE_PRCOMM_DOUBLEV_TAG2,
    EXCHANGE_PRCOMM_NPACKV_TAG,
    EXCHANGE_PRCOMM_NUNPACKV_TAG,
    UPDATE_ADDRESS_TAG
};

int main()
{
    std::cout << REPLACE_DATA << std::endl;
    std::cout << ADD_DATA << std::endl;
    std::cout << MIN_DATA << std::endl;
    std::cout << MAX_DATA << std::endl;
    std::cout << REPLACE_TRANSLATE_DATA << std::endl;
    std::cout << REPLACE_ROTATE_DATA << std::endl;
    std::cout << ADD_ROTATE_DATA << std::endl;
    std::cout << BITWISE_OR_DATA << std::endl;
    std::cout << MIN_NO_PERIODIC_DATA << std::endl;
    std::cout << MAX_NO_PERIODIC_DATA << std::endl;
    std::cout << ADD_NO_PERIODIC_DATA << std::endl;
    std::cout << ADD_TRANSLATE_DATA << std::endl;
    std::cout << SUBTRACT_DATA << std::endl;
    std::cout << SUBTRACT_ROTATE_DATA << std::endl;
    std::cout << NO_CHANGE_DATA << std::endl;
    std::cout << UPDATE_R1_TAG << std::endl;
    std::cout << UPDATE_R2_TAG  << std::endl;
    std::cout << UPDATE_R2_REVERSE_TAG  << std::endl;
    std::cout << EXCHANGE_INT_TAG  << std::endl;
    std::cout << EXCHANGE_INT_REVERSE_TAG  << std::endl;
    std::cout << UPDATE_I1_TAG << std::endl;
    std::cout << UPDATE_I1_REVERSE_TAG  << std::endl;
    std::cout << UPDATE_I2_TAG << std::endl;
    std::cout << UPDATE_I2_REVERSE_TAG << std::endl;
    std::cout << UPDATE_R3_TAG << std::endl;
    std::cout << UPDATE_SYMMETRIC_R3_TAG << std::endl;
    std::cout << EXCHANGE_PRCOMM_INTV_TAG1 << std::endl;
    std::cout << EXCHANGE_PRCOMM_INTV_TAG2 << std::endl;
    std::cout << EXCHANGE_PRCOMM_DOUBLEV_TAG1 << std::endl;
    std::cout << EXCHANGE_PRCOMM_DOUBLEV_TAG2 << std::endl;
    std::cout << EXCHANGE_PRCOMM_NPACKV_TAG << std::endl;
    std::cout << EXCHANGE_PRCOMM_NUNPACKV_TAG << std::endl;
    std::cout << UPDATE_ADDRESS_TAG << std::endl;
}
