package LULESH.Comm is

   --- // lulesh-comm

   --x void CommRecv(Domain& domain, Int_t msgType, Index_t xferFields,
   --x               Index_t dx, Index_t dy, Index_t dz,
   --x               bool doRecv, bool planeOnly);
   procedure CommRecv
     (domainn    : not null access Domain_Record;
      msgType    : in Int_t;
      xferFields : in Index_Type;
      dx         : in Index_Type;
      dy         : in Index_Type;
      dz         : in Index_Type;
      doRecv     : in Boolean;
      planeOnly  : in Boolean);

   --x void CommSend(Domain& domain, Int_t msgType,
   --x               Index_t xferFields, Domain_member *fieldData,
   --x               Index_t dx, Index_t dy, Index_t dz,
   --x               bool doSend, bool planeOnly);
   procedure CommSend
     (domainn    : not null access Domain_Record;
      msgType    : in Int_t;
      xferFields : in Index_Type;
      fieldData  : in Domain_member;
      dx         : in Index_Type;
      dy         : in Index_Type;
      dz         : in Index_Type;
      doSend     : in Boolean;
      planeOnly  : in Boolean);

   --x void CommSBN(Domain& domain, Int_t xferFields, Domain_member *fieldData);
   procedure CommSBN
     (domainn    : not null access Domain_Record;
      xferFields : in Int_t;
      fieldData  : out Domain_member);

   --x void CommSyncPosVel(Domain& domain);
   procedure CommSyncPosVel
     (domainn : not null access Domain_Record);

   --x void CommMonoQ(Domain& domain);
   procedure CommMonoQ
     (domainn : not null access Domain_Record);


end LULESH.Comm;
