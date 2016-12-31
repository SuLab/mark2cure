SELECT  ANY_VALUE(`combined_documents`.`document_pk`) as `document_pk`,
        `document_document`.`document_id` as `document_id`,
        `document_document`.`title` as `title`,

        MAX(`combined_documents`.`total_document_relationships`) as `total_document_relationships`,
        MIN(`combined_documents`.`user_document_relationships`) as `user_document_relationships`,

        MAX(`combined_documents`.`community_answered`) as `community_answered`,
        MAX(`combined_documents`.`community_completed`) as `community_completed`,
        MAX(`combined_documents`.`community_progress`) as `community_progress`,

        MAX(`combined_documents`.`user_completed`) as `user_completed`,
        MAX(`combined_documents`.`user_progress`) as `user_progress`,
        MAX(`combined_documents`.`user_answered`) as `user_answered`,
        MAX(`combined_documents`.`user_view_completed`) as `user_view_completed`

FROM (

  /* This 1st table response are all RelationGroup documents that have at
     least 1 relationship annotation (exclusive JOINs on these)

     (Step 2) Group (document_pk, relation_id, relation_unique_contributors, user_contributed)
     on document_pk so we can build full PMID understanding of community and user specific
     progress */
  SELECT  `computed_attempts`.`document_pk` as `document_pk`,
          `computed_attempts`.`total_document_relationships` as `total_document_relationships`,

          /*  This might have < 0.0 but I don't think we use it anywhere (they get filtered out with
              the user_completed boolean) and I don't honestly want to put this whole thing in an
              IF over limit statement

              (TODO) Compare if the relationships the community has completed overlap with the relatoinships
                     the user has completed
          */
          IF( `computed_attempts`.`total_document_relationships` > @user_work_max,
            /* If there were more than 20 relationships, the user remaining is on those 20 */
            (@user_work_max - `computed_attempts`.`relationships_community_completed`) - `computed_attempts`.`total_user_completed`,

            /* If there were 20 or less relationships, the user remaining is on them all */
            (`computed_attempts`.`total_document_relationships` - `computed_attempts`.`relationships_community_completed`) - `computed_attempts`.`total_user_completed`
          ) as `user_document_relationships`,

          /* Community scope */
          `computed_attempts`.`community_answered` as `community_answered`,
          IF(
            `computed_attempts`.`community_answered` = `computed_attempts`.`total_document_relationships`*@k_max
          , TRUE, FALSE) as `community_completed`,
          /*  This is protected from >1 b/c of the K limit per unique
              contributer limit we in the subquery */
          `computed_attempts`.`community_answered`/(`computed_attempts`.`total_document_relationships`*@k_max) as `community_progress`,

          /* User specific scope */
          IF( `computed_attempts`.`total_document_relationships` > @user_work_max,
            /* If there were more than 20 relationships, the user only had to do 20 */
            IF(`computed_attempts`.`total_user_completed` >= @user_work_max, TRUE, FALSE),

            /* If there were 20 or less relationships, the user only had to do them all */
            IF(`computed_attempts`.`total_user_completed` = `computed_attempts`.`total_document_relationships`, TRUE, FALSE)
          ) as `user_completed`,

          /*  This might have > 1.0 but I don't think we use it anywhere and
              I don't honestly want to put this whole thing in an IF over limit statement */
          IF( `computed_attempts`.`total_document_relationships` > @user_work_max,
            /* If there were more than 20 relationships, the user progress is on those 20 */
            `computed_attempts`.`total_user_completed`/@user_work_max,

            /* If there were 20 or less relationships, the user progress is on them all */
            `computed_attempts`.`total_user_completed`/`computed_attempts`.`total_document_relationships`
          ) as `user_progress`,

          `computed_attempts`.`total_user_completed` as `user_answered`,
          `computed_attempts`.`user_view_completed` as `user_view_completed`

  FROM (
    SELECT  `document_relationships`.`document_pk` as `document_pk`,

            (
              SELECT COUNT(`relation_relation`.`id`)
              FROM `relation_relation`
              WHERE `relation_relation`.`document_id` = `document_relationships`.`document_pk`
            ) as `total_document_relationships`,

            SUM(`document_relationships`.`user_contributed`) as `total_user_completed`,

            /* Community scope */
            SUM(`document_relationships`.`relationship_community_completed`) as `relationships_community_completed`,
            SUM(`document_relationships`.`relation_unique_contributors`) as `community_answered`,

            COALESCE((
              /*  Filtering out uncompleted relationship work
                  that was flagged as being a completed view because
                  the logic (prevous + existing) can be messy and
                  we NEVER want a user to be rewarded a "Relationship
                  Document Quest" bonus more than once

                  I coalesce on a subquery with the `document_view`.`completed`
                  b/c it's cleaner to me than joining and filtering views that
                  may not exists (NULL) and uncompleted views */
                SELECT `document_view`.`completed`
                FROM `document_view`
                INNER JOIN `document_section`
                    ON `document_section`.`id` = `document_view`.`section_id`
                WHERE `document_section`.`document_id` = `document_relationships`.`document_pk`
                    AND `document_view`.`task_type` = 'ri'
                    AND `document_view`.`user_id` = @user_id
                    AND `document_view`.`completed` = True

                /*  There may be multiple Views, order first
                    by the completed ones so the LIMIT section
                    is on that */
                ORDER BY `document_view`.`completed` DESC
                LIMIT 1
            ), FALSE) as `user_view_completed`

    FROM (
      /* (Step 1) Select all Documents with any RelationAnnotation instances
         associated with them. Group by Relations so we get aggregated
         user information for specific comparisons. Results have duplicated
         Document PKs as we will add up the child relationships in Step 2
      */

      SELECT  `document_document`.`id` as `document_pk`,
              `relation_relation`.`id` as `relation_id`,

              /*  Counting the unique users included on a relationship level */
              CASE
                  WHEN COUNT(DISTINCT `document_view`.`user_id`) > @k_max THEN @k_max
                  ELSE COUNT(DISTINCT `document_view`.`user_id`)
              END as `relation_unique_contributors`,

              CASE
                  WHEN COUNT(DISTINCT `document_view`.`user_id`) >= @k_max THEN TRUE
                  ELSE FALSE
              END as `relationship_community_completed`,

              /*  If the provided users of interest has provided
                  an answer for the specific relationship */
              IF(
                SUM(`document_view`.`user_id` = @user_id)
              , TRUE, FALSE) as `user_contributed`

      FROM `document_document`

      INNER JOIN `relation_relation`
          ON `relation_relation`.`document_id` = `document_document`.`id`

      INNER JOIN `relation_relationannotation`
          ON `relation_relationannotation`.`relation_id` = `relation_relation`.`id`

      INNER JOIN `document_annotation`
          ON (`document_annotation`.`object_id` = `relation_relationannotation`.`id`
            AND `document_annotation`.`content_type_id` = @rel_ann_content_type_id)

      INNER JOIN `document_view`
          ON `document_view`.`id` = `document_annotation`.`view_id`

      GROUP BY `relation_id`
    ) as `document_relationships`

    GROUP BY `document_relationships`.`document_pk`

  ) as `computed_attempts`

  UNION ALL

  /* This 2nd Table is all of the RelationGroup documents with blank data
     records with actual annotation data will be overwritten by the first
     table in the union, but allows virgin documents to be included */
  SELECT  `relationship_document`.`id` as `document_pk`,
          @total_available_relations:=(
            SELECT COUNT(`relation_relation`.`id`)
            FROM `relation_relation`
            WHERE `relation_relation`.`document_id` = `relationship_document`.`id`
          )  as `total_document_relationships`,

          IF( @total_available_relations > @user_work_max,
            /* If there were more than 20 relationships, the user remaining is on those 20 */
            @user_work_max,
            /* If there were 20 or less relationships, the user remaining is on them all */
            @total_available_relations
          ) as `user_document_relationships`,

          0 as `community_answered`,
          0 as `community_completed`,
          0 as `community_progress`,

          0 as `user_completed`,
          0 as `user_progress`,
          0 as `user_answered`,
          0 as `user_view_completed`

  FROM (
        /*  Return back the PK, PMID and Title for all Documents
            that are in active RelationGroups */
        SELECT  `document_document`.`id`

        FROM `document_document`

        INNER JOIN `relation_relationgroup_documents`
            ON `relation_relationgroup_documents`.`document_id` = `document_document`.`id`

        INNER JOIN `relation_relationgroup`
            ON `relation_relationgroup`.`id` = `relation_relationgroup_documents`.`relationgroup_id`

        WHERE `relation_relationgroup`.`enabled` = 1
  ) as `relationship_document`

) as `combined_documents`

/* Required to get the PMID + Title at the very end */
LEFT JOIN `document_document`
  ON `document_document`.`id` = `combined_documents`.`document_pk`

/* This doesn't explicitly state to take the first table from the union's results
   but we're using MAX() to get the highest value */
GROUP BY `combined_documents`.`document_pk`

HAVING `community_completed` = FALSE
    AND NOT `user_completed` = TRUE
    AND NOT `user_view_completed` = TRUE
    AND NOT `total_document_relationships` = 0

ORDER BY  `community_progress` DESC,
          `user_document_relationships` DESC

LIMIT 25;
