
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//! Backward search for a character c in an \f$\omega\f$-interval \f$[\ell..r]\f$ in the CSA.
/*!
 * \tparam t_csa CSA type.
 *
 * \param csa    The CSA object.
 * \param l      Left border of the interval \f$ [\ell..r]\f$.
 * \param r      Right border of the interval \f$ [\ell..r]\f$.
 * \param c      Character to be prepended to \f$\omega\f$.
 * \param l_res  New left border.
 * \param r_res  Right border.
 * \return The size of the new interval [\ell_{new}..r_{new}].
 *         Equals zero, if no match is found.
 *
 * \pre \f$ 0 \leq \ell \leq r < csa.size() \f$
 *
 * \par Time complexity
 *         \f$ \Order{ t_{rank\_bwt} } \f$
 * \par Reference
 *         Paolo Ferragina, Giovanni Manzini:
 *         Opportunistic Data Structures with Applications.
 *         FOCS 2000: 390-398
 */
template<class t_csa>
typename t_csa::size_type b_search(
    const t_csa& csa,
    typename t_csa::size_type l,
    typename t_csa::size_type r,
    typename t_csa::char_type c,
    typename t_csa::size_type& l_res,
    typename t_csa::size_type& r_res,
    SDSL_UNUSED typename std::enable_if<std::is_same<csa_tag, typename t_csa::index_category>::value, csa_tag>::type x = csa_tag()
)
{
    assert(l <= r); assert(r < csa.size());
    typename t_csa::size_type cc = csa.char2comp[c];
    if (cc == 0 and c > 0) {
        l_res = 1;
        r_res = 0;
    } else {
        typename t_csa::size_type c_begin = csa.C[cc];
        if (l == 0 and r + 1 == csa.size()) {
            l_res = c_begin;
            r_res = csa.C[cc + 1] - 1;
        } else {
            l_res = c_begin + csa.bwt.rank(l, c); // count c in bwt[0..l-1]
            r_res = c_begin + csa.bwt.rank(r + 1, c) - 1; // count c in bwt[0..r]
        }
    }

    ///My printing: should be removed - enam

    // std:: cout << "saa - backward_search : " << c << " l = " << l << ", r = " << r << ", l_res = " << l_res << ", r_res = " << r_res << ", final = " << (r_res + 1 - l_res) << std:: endl;

    assert(r_res + 1 - l_res >= 0);
    return r_res + 1 - l_res;
}
//! Backward search for a pattern in an \f$\omega\f$-interval \f$[\ell..r]\f$ in the CSA.
/*!
 * \tparam t_csa      A CSA type.
 * \tparam t_pat_iter Pattern iterator type.
 *
 * \param csa   The CSA object.
 * \param l     Left border of the lcp-interval \f$ [\ell..r]\f$.
 * \param r     Right border of the lcp-interval \f$ [\ell..r]\f$.
 * \param begin Iterator to the begin of the pattern (inclusive).
 * \param end   Iterator to the end of the pattern (exclusive).
 * \param l_res New left border.
 * \param r_res New right border.
 * \return The size of the new interval [\ell_{new}..r_{new}].
 *         Equals zero, if no match is found.
 *
 * \pre \f$ 0 \leq \ell \leq r < csa.size() \f$
 *
 * \par Time complexity
 *       \f$ \Order{ len \cdot t_{rank\_bwt} } \f$
 * \par Reference
 *         Paolo Ferragina, Giovanni Manzini:
 *         Opportunistic Data Structures with Applications.
 *         FOCS 2000: 390-398
 */
template<class t_csa, class t_pat_iter>
typename t_csa::size_type
b_search(
    const t_csa& csa,
    typename t_csa::size_type l,
    typename t_csa::size_type r,
    t_pat_iter begin,
    t_pat_iter end,
    typename t_csa::size_type& l_res,
    typename t_csa::size_type& r_res,
    typename t_csa::size_type& min_k,
    typename t_csa::size_type& max_count,
    typename t_csa::size_type& k_length,
    SDSL_UNUSED typename std::enable_if<std::is_same<csa_tag, typename t_csa::index_category>::value, csa_tag>::type x = csa_tag()
)
{
    t_pat_iter it = end;
    /// saves the last results
    k_length = 0;
    typename t_csa::size_type last_l = -1, last_r = -1,  res;
    // typename t_csa::size_type last_l_res = -1, last_r_res = -1;
    while (begin < it and r + 1 - l > 0) {
        last_l = l;
        last_r = r;
        --it;
        b_search(csa, l, r, (typename t_csa::char_type)*it, l, r);
        // last_l_res = l;
        // last_r_res = r;
        k_length++;
    }
    /// ended at the beginning
    if (r + 1 - l > 0)
    {
        l_res = l;
        r_res = r;
        res = r + 1 - l;
        if (res <= max_count && k_length >= min_k) return res;
        return 0;
    }
    k_length--;
    if (k_length < min_k) return 0;

    // cout << "k_length = " << k_length << endl;
    l_res = last_l;
    r_res = last_r;
    res = last_r + 1 - last_l;
    if (res <= max_count && k_length >= min_k) return res;
    return 0;
}
//! Calculates all occurrences of a pattern pat in a CSA.
/*!
 * \tparam t_csa      CSA type.
 * \tparam t_pat_iter Pattern iterator type.
 * \tparam t_rac      Resizeable random access container.
 *
 * \param csa   The CSA object.
 * \param begin Iterator to the begin of the pattern (inclusive).
 * \param end   Iterator to the end of the pattern (exclusive).
 * \return A vector containing the occurrences of the pattern  in the CSA.
 *
 * \par Time complexity
 *        \f$ \Order{ t_{backward\_search} + z \cdot t_{SA} } \f$, where \f$z\f$ is the number of
 *         occurrences of pattern in the CSA.
 */
template<class t_csa, class t_pat_iter, class t_rac = int_vector<64>>
t_rac my_locate(
    const t_csa&  csa,
    t_pat_iter begin,
    t_pat_iter end,
    typename t_csa::size_type& min_k,
    typename t_csa::size_type& max_count,
    typename t_csa::size_type& k_length,
    SDSL_UNUSED typename std::enable_if<std::is_same<csa_tag, typename t_csa::index_category>::value, csa_tag>::type x = csa_tag()
)
{
    // cout << "min(k) = " <<  min_k << " max_count = " << max_count << endl;

    typename t_csa::size_type occ_begin, occ_end, occs;

    occs = b_search(csa, 0, csa.size() - 1, begin, end, occ_begin, occ_end, min_k, max_count, k_length);

    /// it should be handled in b_search() function. Kept it here for safety. it should be removed later.
    if ( (k_length < min_k && occs > 0) || occs > max_count)
    {
        cerr << "It should not come here!" << endl;

        t_rac occ(0);
        return occ;
    }
    t_rac occ(occs);
    for (typename t_csa::size_type i = 0; i < occs; ++i) {
        occ[i] = csa[occ_begin + i];
    }
    return occ;
}
//! Calculates all occurrences of a pattern pat in a CSA/CST.
/*!
 * \tparam t_csa      CSA/CST type.
 * \tparam t_rac      Resizeable random access container.
 *
 * \param csa  The CSA/CST object.
 * \param pat  The pattern.
 * \return A vector containing the occurrences of the pattern  in the CSA.
 *
 * \par Time complexity
 *        \f$ \Order{ t_{backward\_search} + z \cdot t_{SA} } \f$, where \f$z\f$ is the number of
 *         occurrences of pattern in the CSA.
 */
// template<class t_csx, class t_rac = int_vector<64>>
// t_rac my_locate(
//     const t_csx&  csx,
//     const typename t_csx::string_type& pat,
//     typename t_csx::size_type& min_k,
//     typename t_csx::size_type& max_count,
//     typename t_csx::size_type& k_length
// )
// {
//     typename t_csx::index_category tag;
//     return my_locate<t_csx, decltype(pat.begin()), t_rac>(csx, pat.begin(), pat.end(), tag, min_k, max_count, k_length);
// }
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>