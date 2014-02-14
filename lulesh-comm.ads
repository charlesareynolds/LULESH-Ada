package LULESH.Comm is

   --- // lulesh-comm

   --x void CommRecv(Domain& domain, Int_t msgType, Index_t xferFields,
   --x               Index_t dx, Index_t dy, Index_t dz,
   --x               bool doRecv, bool planeOnly);
   procedure CommRecv
     (domainn    : not null access Domain_Record;
      msgType    : in Int_t;
      xferFields : in Index_t;
      dx         : in Index_t;
      dy         : in Index_t;
      dz         : in Index_t;
      doRecv     : in Boolean;
      planeOnly  : in Boolean);

   --x void CommSend(Domain& domain, Int_t msgType,
   --x               Index_t xferFields, Domain_member *fieldData,
   --x               Index_t dx, Index_t dy, Index_t dz,
   --x               bool doSend, bool planeOnly);
   procedure CommSend
     (domainn    : not null access Domain_Record;
      msgType    : in IC.int;
      xferFields : in Index_t;
      fieldData  : in Domain_member;
      dx         : in Index_t;
      dy         : in Index_t;
      dz         : in Index_t;
      doSend     : in Boolean;
      planeOnly  : in Boolean);

   --x void CommSBN(Domain& domain, Int_t xferFields, Domain_member *fieldData);
   procedure CommSBN
     (domainn    : not null access Domain_Record;
      xferFields : in Index_t;
      fieldData  : out Domain_member);

   --x void CommSyncPosVel(Domain& domain);
   procedure CommSyncPosVel
     (domainn : not null access Domain_Record);

   --x void CommMonoQ(Domain& domain);
   procedure CommMonoQ
     (domainn : not null access Domain_Record);


end LULESH.Comm;
