C/*****************************************************************************/
C/*                                                                           */
C/*   (C) Copyright IBM Corporation, 1997                                     */
C/*   (C) Copyright Regents of the University of Minnesota, 1997              */
C/*                                                                           */
C/*   genleaftree.f                                                           */
C/*                                                                           */
C/*   Written by Mahesh Joshi, U of MN.                                       */
C/*                                                                           */
C/*****************************************************************************/
C/*                                                                           */
C/* This code is meant to be used solely for educational, research, and       */
C/* benchmarking purposes by non-profit institutions and US government        */
C/* agencies only.  Use by any other organization requires prior written      */
C/* permission from both IBM Corporation and the University of Minnesota.     */
C/* The software may not be sold or redistributed.  One may make copies       */
C/* of the software or modify it for their use provided that the copies,      */
C/* modified or otherwise, are not sold or distributed, are used under the    */
C/* same terms and conditions, and this notice and any part of the source     */
C/* code that follows this notice are not separated.                          */
C/*                                                                           */
C/* As unestablished research software, this code is provided on an           */
C/* ``as is'' basis without warranty of any kind, either expressed or         */
C/* implied, including but not limited to implied warranties of               */
C/* merchantability and fitness for a particular purpose.  IBM does not       */
C/* warrant that the functions contained in this software will meet the       */
C/* user's requirements or that the operation of its routines will be         */
C/* uninterrupted or error-free.  Acceptance and use of this program          */
C/* constitutes the user's understanding that he/she will have no recourse    */
C/* to IBM for any actual or consequential damages, including, but not        */
C/* limited to, lost profits or savings, arising out of the use or inability  */
C/* to use these libraries.  Even if the user informs IBM of the possibility  */
C/* of such damages, IBM expects the user to accept the risk of any such      */
C/* harm, or the user shall not attempt to use these libraries for any        */
C/* purpose.                                                                  */
C/*                                                                           */
C/* The downloading, compiling, or executing any part of this software        */
C/* constitutes an implicit agreement to these terms.  These terms and        */
C/* conditions are subject to change at any time without prior notice.        */
C/*                                                                           */
C/*****************************************************************************/
C/* $Id: genleaftree.f,v 1.1 2004-12-10 20:26:45 paklein Exp $ */
C/*****************************************************************************/

      subroutine GENLEAFTREE(beginnode,size,parent,ptrs,inds,setid)

      integer beginnode, size
      integer parent(0:*), ptrs(2,0:*), inds(*), setid(0:*)
      integer i,j,k,l,lset,node,endnode,tpar

      tpar = -1

      do 30 i = beginnode, beginnode+size-1
        parent(i) = tpar
        setid(i) = tpar
        k = ptrs(1,i)
        do 40 j = k, k + ptrs(2,i) - 1
          l = inds(j)

          if (l .ge. i) go to 30

          lset = l
          node = setid(l)
          do while (node .ne. tpar)
            setid(lset) = i
            lset = node
            node = setid(node)
          end do
      
          if(lset.ne.i) then
            parent(lset) = i
            setid(lset) = i
          end if

40      continue
30    continue

      end
