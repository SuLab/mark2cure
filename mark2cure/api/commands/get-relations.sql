SELECT  `document_relationships`.`document_pk`,
        `document_relationships`.`pmid`,
        `document_relationships`.`title`,

    COUNT(DISTINCT `document_relationships`.`relation_id`) as `document_relationships`,
    SUM(`document_relationships`.`community_completed`)/COUNT(DISTINCT `document_relationships`.`relation_id`) as `community_progress`,

    /* (TODO) The Community and User Completed boolean
        can't be based off if their relationship submissions
        is equal to that in the docuemnt, but should be based
        of a static number.

        This supports the feature to prevent having to limit
        documents b/c there are too many relationship
        comparisons within it */

    /*  Yes it's redundant with the progress but the
        boolean is nice to use */
    IF(
      SUM(`document_relationships`.`community_completed`) = COUNT(DISTINCT `document_relationships`.`relation_id`)
    , TRUE, FALSE) as `community_completed`,

    IF(
      SUM(`document_relationships`.`user_completed`) = COUNT(DISTINCT `document_relationships`.`relation_id`)
    , TRUE, FALSE) as `user_completed`

FROM (
  SELECT  `relationship_document`.`id` as `document_pk`,
          `relationship_document`.`document_id` as `pmid`,
          `relationship_document`.`title` as `title`,

          `relation_relation`.`id` as `relation_id`,

          # (TODO) What happends when there are no annotations?
          # MIN(`document_annotation`.`created`) as `first_annotation`,
          # MAX(`document_annotation`.`created`) as `last_annotation`,

      /*  If the relationship has been answered enough
          times by a minimum number of unique users */
      IF(
        /*  Ensure we count the unique number of users who
            have completed the relationship */
        COUNT(DISTINCT `document_view`.`user_id`) >= {completions}
          , TRUE, FALSE) as `community_completed`,

      /*  If the provided users of interest has provided
          an answer for the specific relationship */
          IF(
            SUM(`document_view`.`user_id` = {user_id})
          , TRUE, FALSE) as `user_completed`

  FROM (
    /* Return back the PK, PMID and Title for all Documents
       that are in active RelationGroups */
      SELECT  `document_document`.`id` as `id`,
              `document_document`.`document_id`,
              `document_document`.`title` as `title`

      FROM `document_document`

      INNER JOIN `relation_relationgroup_documents`
          ON `relation_relationgroup_documents`.`document_id` = `document_document`.`id`

      INNER JOIN `relation_relationgroup` as `group`
          ON `group`.`id` = `relation_relationgroup_documents`.`relationgroup_id`

      WHERE `group`.`enabled` = 1

  ) as `relationship_document`

  /*  This INNER join is somewhat redundant as enforces only
      Documents with known potential relationships to show up.
      However, I'm still keeping the subquery as I think there is value
      to selectively fetch documents from a RelationGroup that have
      it's own respective controls (eg. Enabled) */
  INNER JOIN `relation_relation`
      ON `relation_relation`.`document_id` = `relationship_document`.`id`

  LEFT JOIN `relation_relationannotation`
      ON `relation_relationannotation`.`relation_id` = `relation_relation`.`id`

  LEFT JOIN `document_annotation`
      ON (`document_annotation`.`object_id` = `relation_relationannotation`.`id`
        AND `document_annotation`.`content_type_id` = 56)

  /*  We're only joining Views to get back the User ID,
      We don't care about the "completed" boolean b/c
      the Relation tasks are not all or none */
  INNER JOIN `document_view`
      ON `document_view`.`id` = `document_annotation`.`view_id`

  GROUP BY `relation_id`
) as `document_relationships`

GROUP BY `document_relationships`.`document_pk`

HAVING `community_completed` = FALSE AND NOT `user_completed` = TRUE

ORDER BY  `community_progress` DESC,
          `document_relationships` DESC

LIMIT 25
